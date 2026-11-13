"""Deja un ligando listo para GROMACS, en la pose que le dio GNINA.

Hace en un paso lo que en la simulacion del 453 fueron cinco scripts. Alli el
SDF tenia el confomero de la biblioteca y el mol2 la pose dockeada, con los
atomos en ordenes distintos, asi que habia que mapear uno contra otro y
reconstruir el .gro a mano. Aqui se construye primero un SDF que YA tiene la
pose dockeada y los hidrogenos, y entonces el .gro que exporta Interchange es
directamente la pose buena.

Lo que hace, en orden:

  1. Lee los ordenes de enlace del representante estructural (el
     registro <N> de representantes_docking.sdf, sin hidrogenos) y las
     coordenadas de la pose 1 del mol2 de GNINA.
  2. Empareja los atomos de los dos ficheros. Vienen en ordenes distintos, y
     las etiquetas canonicas no bastan porque las moleculas con simetria tienen
     atomos indistinguibles; se resuelve con un isomorfismo de grafos.
  3. Anade los hidrogenos por geometria y los relaja con MMFF manteniendo fijos
     los atomos pesados, que son los que trae el docking.
  4. Parametriza con OpenFF 2.2.1 y cargas NAGL, el mismo campo de fuerza que
     se uso con el 453, y exporta a GROMACS con topologia modular.
  5. Comprueba que el .gro exportado conserva la pose: el RMSD de los atomos
     pesados contra el mol2 tiene que salir 0,000 A.

El ligando se pasa como argumento (por defecto, 170):

    conda run -n openff python scripts/md/ligando/preparar_ligando.py 170

Entorno: `openff` (openff-toolkit, openff-interchange, openff-nagl, rdkit).

Entrada, relativa a la raiz del repositorio:
  results/docking/representantes_docking.sdf   (registro <N>: los ordenes de
      enlace; el numero del ligando es su posicion en ese fichero)
  results/docking/candidatos/docked_ligando_<N>.mol2   (la pose de GNINA)
Salida: WORK/ligand/ (por defecto run/ligand/), y la copia versionada de la
ejecucion original esta en data/md/ligand_<N>/gromacs/.
"""

import os
import re
import sys
from pathlib import Path

import numpy as np
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule
from openff.units import unit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

LIGANDO = sys.argv[1] if len(sys.argv) > 1 else "170"
RAIZ = Path(os.environ.get("RAIZ", "."))
WORK = Path(os.environ.get("WORK", "run"))
ENTRADA = Path(os.environ.get("ENTRADA", RAIZ / "results" / "docking"))
SALIDA = WORK / "ligand"
REPRESENTANTES = ENTRADA / "representantes_docking.sdf"
MOL2 = ENTRADA / "candidatos" / f"docked_ligando_{LIGANDO}.mol2"
POSE = 1                       # la pose que selecciono el cribado
PREFIJO = f"ligando_{LIGANDO}"
NAGL = "openff-gnn-am1bcc-1.0.0.pt"
CAMPO = "openff-2.2.1.offxml"


# --- lectura ------------------------------------------------------------------
def lee_mol2(path, pose=1):
    """Una pose del mol2 de GNINA.

    Entre ATOM y BOND puede haber secciones intercaladas (algunos ficheros
    traen UNITY_ATOM_ATTR), asi que cada seccion se corta en el siguiente
    @<TRIPOS> y las lineas se filtran por forma.
    """
    t = path.read_text().split("@<TRIPOS>MOLECULE")[pose]

    def seccion(nombre):
        return [l.split() for l in t.split(f"@<TRIPOS>{nombre}")[1]
                .split("@<TRIPOS>")[0].strip().splitlines() if l.strip()]

    at = [f for f in seccion("ATOM") if len(f) >= 6 and f[0].isdigit()]
    bo = [f for f in seccion("BOND") if len(f) >= 4 and f[0].isdigit()]
    els = [re.sub(r"[^A-Za-z]", "", f[5].split(".")[0]) for f in at]
    X = np.array([[float(x) for x in f[2:5]] for f in at])
    bonds = [(int(f[1]) - 1, int(f[2]) - 1) for f in bo]
    return els, X, bonds


def grafo_solo_conectividad(els, bonds):
    """Molecula con todos los enlaces sencillos y sin sanear, para que el
    emparejamiento de subestructura sea un isomorfismo de grafos sobre
    elementos y no estorben aromaticidad ni valencias."""
    m = Chem.RWMol()
    for e in els:
        m.AddAtom(Chem.Atom(e))
    for a, b in bonds:
        m.AddBond(a, b, Chem.BondType.SINGLE)
    mm = m.GetMol()
    Chem.SanitizeMol(mm, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
    Chem.FastFindRings(mm)
    return mm


print("=" * 62)
print(f"PREPARACION DEL LIGANDO  ({PREFIJO}, pose {POSE} de GNINA)")
print("=" * 62)

def representante(n):
    """El registro n-esimo (1-indexado) de representantes_docking.sdf.

    El numero del ligando ES su posicion en ese fichero, que es el orden en
    que lo dejo el cribado: el 170 es el registro 170. Los .mol2 de
    candidatos/ llevan ese mismo numero en el nombre.
    """
    registros = [r for r in REPRESENTANTES.read_text().split("$$$$") if r.strip()]
    if not 1 <= n <= len(registros):
        sys.exit(f"{REPRESENTANTES} tiene {len(registros)} registros; pedido el {n}")
    return registros[n - 1].lstrip("\n")


bloque = representante(int(LIGANDO))
ref = Chem.MolFromMolBlock(bloque, removeHs=False)
if ref is None:
    sys.exit(f"RDKit no puede leer el registro {LIGANDO} de {REPRESENTANTES}")
nombre = bloque.splitlines()[0].split()[0]
print(f"\n1. Representante estructural: {nombre}")
print(f"   {ref.GetNumAtoms()} atomos pesados, "
      f"{Chem.rdMolDescriptors.CalcMolFormula(ref)}")
print(f"   SMILES: {Chem.MolToSmiles(ref)}")

els_m, X_m, bonds_m = lee_mol2(MOL2, POSE)
pes_m = [i for i, e in enumerate(els_m) if e != "H"]
print(f"\n2. Pose de GNINA: {len(els_m)} atomos, {len(pes_m)} pesados")

# --- emparejamiento -----------------------------------------------------------
im = {a: k for k, a in enumerate(pes_m)}
gA = grafo_solo_conectividad(
    [a.GetSymbol() for a in ref.GetAtoms()],
    [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in ref.GetBonds()])
gB = grafo_solo_conectividad(
    [els_m[a] for a in pes_m],
    [(im[a], im[b]) for a, b in bonds_m
     if els_m[a] != "H" and els_m[b] != "H"])
match = gB.GetSubstructMatch(gA)
if len(match) != ref.GetNumAtoms():
    sys.exit(f"isomorfismo fallido: {len(match)}/{ref.GetNumAtoms()} atomos")
print(f"   emparejados {len(match)}/{ref.GetNumAtoms()} atomos pesados")

conf = ref.GetConformer()
for k in range(ref.GetNumAtoms()):
    conf.SetAtomPosition(k, Point3D(*X_m[pes_m[match[k]]]))

# --- hidrogenos ---------------------------------------------------------------
molH = Chem.AddHs(ref, addCoords=True)
mp = AllChem.MMFFGetMoleculeProperties(molH)
if mp is None:
    sys.exit("MMFF no sabe tipar esta molecula; revisa el SDF de entrada")
ff_mmff = AllChem.MMFFGetMoleculeForceField(molH, mp)
for i in range(ref.GetNumAtoms()):
    ff_mmff.AddFixedPoint(i)          # los pesados no se tocan: son el docking
antes = ff_mmff.CalcEnergy()
ff_mmff.Minimize(maxIts=2000)
print(f"\n3. Hidrogenos: {molH.GetNumAtoms() - ref.GetNumAtoms()} anadidos, "
      f"relajados con los pesados fijos "
      f"(MMFF {antes:.1f} -> {ff_mmff.CalcEnergy():.1f} kcal/mol)")

SALIDA.mkdir(parents=True, exist_ok=True)
sdf_pose = SALIDA / f"{PREFIJO}.sdf"
w = Chem.SDWriter(str(sdf_pose))
molH.SetProp("_Name", nombre)
w.write(molH)
w.close()
print(f"   escrito {sdf_pose}")

# --- parametrizacion ----------------------------------------------------------
print(f"\n4. Parametrizando con {CAMPO} y cargas NAGL")
lig = Molecule.from_file(str(sdf_pose), file_format="SDF")
# OpenFF usa el nombre de la molecula en la exportacion a GROMACS, y el del
# SDF trae varios identificadores separados por espacios: no vale.
lig.name = "LIG"
lig.assign_partial_charges(partial_charge_method=NAGL,
                           normalize_partial_charges=True)
print(f"   atomos: {lig.n_atoms}   carga total: {lig.total_charge}")

inter = Interchange.from_smirnoff(force_field=ForceField(CAMPO),
                                 topology=lig.to_topology(),
                                 charge_from_molecules=[lig])
# Caja tecnica, solo para poder exportar. La de verdad se define despues
# alrededor del complejo con `gmx editconf`.
inter.box = unit.Quantity(np.diag([2.5, 2.5, 2.5]), unit.nanometer)
inter.to_gromacs(prefix=str(SALIDA / PREFIJO), monolithic=False)
print(f"   exportado a {SALIDA}/{PREFIJO}*")

# --- comprobacion -------------------------------------------------------------
# El .gro tiene que ser la pose del docking, no otra conformacion. Si el RMSD
# de los pesados no sale 0,000 el emparejamiento esta mal y no hay que seguir.
gro = (SALIDA / f"{PREFIJO}.gro").read_text().splitlines()
n = int(gro[1])
X_gro = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                  for l in gro[2:2 + n]]) * 10.0        # nm -> A
els_gro = [a.symbol for a in lig.atoms]
pes_gro = [i for i, e in enumerate(els_gro) if e != "H"]
ref_xyz = np.array([X_m[pes_m[match[k]]] for k in range(ref.GetNumAtoms())])
rmsd = float(np.sqrt(((X_gro[pes_gro] - ref_xyz) ** 2).sum(1).mean()))
print(f"\n5. RMSD de los atomos pesados, .gro contra la pose de GNINA: "
      f"{rmsd:.3f} A")
if rmsd > 0.01:
    sys.exit("El .gro NO conserva la pose dockeada. No sigas.")
print("   correcto: el .gro conserva la pose del docking")

print("\n" + "=" * 62)
print("LIGANDO LISTO")
print("=" * 62)
for f in sorted(SALIDA.iterdir()):
    print(f"  {f}  ({f.stat().st_size} bytes)")
