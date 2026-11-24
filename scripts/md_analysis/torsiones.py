# Barrido de TODOS los enlaces rotables del ligando a lo largo de la
# trayectoria, con deteccion de cambios de rotamero sostenidos.
#
# El objetivo no es medir el enlazador del benzotriazol: es medirlo todo, para
# poder decir que de los diez grados de libertad torsionales del ligando
# exactamente uno cambia en la ventana de la transicion. Medir solo el que uno
# sospecha contesta una pregunta mas debil.
#
# Resultado: la torsion CH2-CH2 del enlazador salta en 27,75 y 30,65 ns. Los
# demas enlaces, o no saltan, o lo hacen a partir de 31,25 ns, cuando el
# reacomodo ya ha terminado. Las dos amidas son los enlaces de menor amplitud.
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/ligand_453/frame0.pdb
#   open data/md/ligand_453/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/torsiones.py
import csv
import os
import sys

import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.geometry import dihedral

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import encuentra_amidas, exige_trayectoria, hidrogenos, vecinos  # noqa: E402

# Ligando y paso entre fotogramas como argumentos:
#   runscript scripts/md_analysis/torsiones.py 170 0.1
# Sin argumentos, el 453 y sus 50 ps. Lee lo que este cargado en #1 y escribe
# en results/md/ligand_<N>/tables/.
LIGANDO = sys.argv[1] if len(sys.argv) > 1 else "453"
DT_NS = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
OUT = f"results/md/ligand_{LIGANDO}/tables/"
os.makedirs(OUT, exist_ok=True)
VENTANA = (27.5, 30.0)          # la transicion, segun results/md/ligand_453/README.md
SOSTENIDO = 20                  # fotogramas (1 ns) a cada lado de un salto
SALTO_MIN = 60                  # grados minimos para contarlo como salto

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
exige_trayectoria(s)   # falla antes de sobrescribir nada
lig = s.residues.filter(s.residues.names == "LIG")[0]

# --- Identificacion quimica de cada posicion ---------------------------------
anillos = []
for ring in s.rings(cross_residue=False, all_size_threshold=6):
    ats = list(ring.ordered_atoms)
    if ats[0].residue is lig:
        anillos.append(ats)
en_anillo = {a for ats in anillos for a in ats}


def nombre_anillo(ats):
    """El anillo, por su composicion. No depende del orden de enumeracion."""
    els = [a.element.name for a in ats]
    if "S" in els:
        return "tiofeno"
    if els.count("N") >= 2:
        return "triazol"
    if len(ats) == 6 and set(els) == {"C"}:
        # El dimetoxifenilo es el unico anillo de seis carbonos con oxigenos
        # colgando; el otro es el benceno fusionado del benzotriazol.
        if any(v.element.name == "O" for a in ats for v in vecinos(a)):
            return "fenilo"
        return "benzo"
    return "anillo"


anillo_de = {a: nombre_anillo(ats) for ats in anillos for a in ats}

# Papel de cada atomo que pertenece a una amida o al enlazador. Es lo que
# permite distinguir dos enlaces que, mirados de cerca, son identicos: el C-N
# de una amida y el C-N de la otra.
papel = {}
amidas = encuentra_amidas(s, lig)
for grupo, c_acilo, c, o, n, c_amino, _ in amidas:
    papel[c] = (grupo, "CO")
    papel[n] = (grupo, "Nam")

# El enlazador: la cadena que sale del carbonilo del benzotriazol hasta el
# anillo, mientras no se bifurque.
bt = [a for a in amidas if a[0] == "benzotriazol"][0]
cadena = [bt[2], bt[1]]                # carbonilo, C_acilo
while True:
    siguientes = [x for x in vecinos(cadena[-1]) if x is not cadena[-2]]
    if len(siguientes) != 1:
        break
    cadena.append(siguientes[0])
for k, a in enumerate(cadena[1:-1]):
    papel[a] = ("enlazador", f"CH{len(hidrogenos(a))}")


def descriptor(a):
    if a in papel:
        return papel[a][1]
    if a in en_anillo:
        return anillo_de[a]
    if a.element.name == "O":
        return "Ometoxi"
    if a.element.name == "C":
        return f"CH{len(hidrogenos(a))}"
    return a.element.name


def grupo_de(a1, a2):
    """El sustituyente al que pertenece el enlace, si se puede decir."""
    for a in (a1, a2):
        if a in papel:
            return papel[a][0]
    return None


# --- Enlaces rotables ---------------------------------------------------------
tors, nombres = [], []
for b in lig.atoms.intra_bonds:
    a1, a2 = b.atoms
    if "H" in (a1.element.name, a2.element.name):
        continue                                   # no es del esqueleto
    if a1 in en_anillo and a2 in en_anillo:
        continue                                   # interno de un anillo
    v1 = [x for x in vecinos(a1) if x is not a2]
    v2 = [x for x in vecinos(a2) if x is not a1]
    if not v1 or not v2:
        continue                                   # un extremo terminal
    # Los sustituyentes de referencia se eligen por indice de atomo: cualquier
    # criterio vale mientras sea reproducible y no dependa del orden en que la
    # herramienta enumere los vecinos.
    ref1 = min(v1, key=lambda a: a.coord_index)
    ref2 = min(v2, key=lambda a: a.coord_index)
    grupo = grupo_de(a1, a2)
    nombre = f"{descriptor(a1)}-{descriptor(a2)}"
    tors.append(((ref1, a1, a2, ref2), grupo, nombre,
                 min(a1.coord_index, a2.coord_index)))

# Nombre final: prefijo de grupo, y un sufijo numerico solo si hace falta para
# desempatar (los dos metoxis son indistinguibles por composicion).
completos = [f"{g + '_' if g else ''}{n}" for _, g, n, _ in tors]
nombres = []
for i, nombre in enumerate(completos):
    repes = [j for j, x in enumerate(completos) if x == nombre]
    nombres.append(nombre if len(repes) == 1 else f"{nombre}{repes.index(i) + 1}")

orden = sorted(range(len(tors)), key=lambda i: (tors[i][1] or "zz", tors[i][3]))
tors = [tors[i][0] for i in orden]
nombres = [nombres[i] for i in orden]
session.logger.info(f"enlaces rotables: {len(tors)}")

# --- Medida -------------------------------------------------------------------
filas = []
for cid in s.coordset_ids:
    s.active_coordset_id = cid
    filas.append([round((cid - 1) * DT_NS, 2)]
                 + [round(dihedral(*(x.coord for x in t)), 1) for t in tors])
s.active_coordset_id = 1

with open(OUT + "torsiones.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["tiempo_ns"] + nombres)
    wr.writerows(filas)

# --- Resumen ------------------------------------------------------------------
t = np.array([f[0] for f in filas])


def media_circular(v):
    return np.degrees(np.arctan2(np.sin(np.radians(v)).mean(),
                                 np.cos(np.radians(v)).mean()))


def diferencia(a, b):
    return abs((b - a + 180) % 360 - 180)


def amplitud(v):
    """Desviacion tipica circular: cuanto se mueve la torsion, sin importar
    alrededor de que valor."""
    ang = np.radians(v)
    r = np.hypot(np.sin(ang).mean(), np.cos(ang).mean())
    return np.degrees(np.sqrt(-2 * np.log(max(r, 1e-12))))


session.logger.info(f"{'enlace':34s} {'amplitud':>9s}  saltos sostenidos (ns)")
for k, nombre in enumerate(nombres):
    x = np.array([f[1 + k] for f in filas])
    d = np.array([diferencia(media_circular(x[i - SOSTENIDO:i]),
                             media_circular(x[i:i + SOSTENIDO]))
                  for i in range(SOSTENIDO, len(x) - SOSTENIDO)])
    tt = t[SOSTENIDO:len(x) - SOSTENIDO]
    # Un salto, no cada fotograma por encima del umbral: solo el maximo local.
    picos = [tt[i] for i in range(len(d))
             if d[i] > SALTO_MIN
             and d[i] == d[max(0, i - SOSTENIDO):i + SOSTENIDO + 1].max()]
    dentro = [p for p in picos if VENTANA[0] <= p <= VENTANA[1]]
    session.logger.info(
        f"{nombre:34s} {amplitud(x):8.1f}°  "
        + (", ".join(f"{p:.2f}" for p in picos) or "—")
        + ("   <<< en la ventana" if dentro else "")
    )
