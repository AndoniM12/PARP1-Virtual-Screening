"""Une proteina y ligando: la geometria y la topologia.

Dos cosas que en la simulacion del 453 se hicieron a mano (un heredoc de Python
para el .gro y un `sed` para la topologia):

  1. Concatena los dos .gro en uno. La caja se toma de la proteina y se
     redefine despues con `gmx editconf`.
  2. Mete el ligando en topol.top: el #include de los tipos de atomo ANTES de
     [ moleculetype ] --ahi es donde GROMACS espera la seccion
     [ atomtypes ]--, el #include del propio ligando despues de los iones, y
     la linea `LIG 1` en [ molecules ].

Es idempotente: si la topologia ya tiene el ligando, no lo mete dos veces.

El ligando se pasa como argumento (por defecto, 170):

    python scripts/md/ligando/montar_complejo.py 170
"""

import os
import sys
from pathlib import Path

LIG_N = sys.argv[1] if len(sys.argv) > 1 else "170"
WORK = Path(os.environ.get("WORK", "run"))
PREFIJO = f"ligando_{LIG_N}"
PROTEINA = WORK / "protein"
LIGANDO = WORK / "ligand"
SALIDA = WORK / "complex" / "complejo.gro"

gro_p = PROTEINA / "PARP1.gro"
gro_l = LIGANDO / f"{PREFIJO}.gro"
top = PROTEINA / "topol.top"

for f in (gro_p, gro_l, top):
    if not f.exists():
        sys.exit(f"falta {f}; ejecuta antes las etapas previas")

# --- 1. geometria -------------------------------------------------------------
lp = gro_p.read_text().splitlines()
ll = gro_l.read_text().splitlines()
np_, nl = int(lp[1].strip()), int(ll[1].strip())
caja = lp[2 + np_]

SALIDA.parent.mkdir(parents=True, exist_ok=True)
with SALIDA.open("w") as f:
    f.write(f"PARP1 - {PREFIJO} complex\n")
    f.write(f"{np_ + nl:5d}\n")
    for l in lp[2:2 + np_]:
        f.write(l + "\n")
    for l in ll[2:2 + nl]:
        f.write(l + "\n")
    f.write(caja + "\n")

print(f"complejo: {SALIDA}")
print(f"  proteina {np_} atomos + ligando {nl} = {np_ + nl}")

# --- 2. topologia -------------------------------------------------------------
# Las rutas van relativas al DIRECTORIO DE topol.top, no al directorio de
# trabajo: el preprocesador de GROMACS busca cada #include junto al fichero que
# lo incluye (y en GMXLIB), no donde se lance grompp. La linea de posre.itp que
# escribe pdb2gmx si lleva la ruta del directorio de trabajo, pero esa va
# dentro de un #ifdef POSRES y no se procesa, asi que no sirve de ejemplo.
rel = Path("..") / "ligand"
inc_tipos = f'#include "{rel / (PREFIJO + "_atomtypes.itp")}"'
inc_lig = f'#include "{rel / (PREFIJO + "_LIG.itp")}"'

lineas = top.read_text().splitlines()
if any(inc_lig == l.strip() for l in lineas):
    print(f"topologia: {top} ya incluye el ligando, sin cambios")
    sys.exit(0)

salida, puesto_tipos, puesto_lig, puesto_mol = [], False, False, False
en_molecules = False
for l in lineas:
    s = l.strip()
    # los tipos de atomo, justo tras el campo de fuerza y antes de moleculetype
    if not puesto_tipos and s.startswith("#include") and "forcefield.itp" in s:
        salida.append(l)
        salida.append(inc_tipos)
        puesto_tipos = True
        continue
    # el ligando, tras el include de los iones
    if not puesto_lig and s.startswith("#include") and "ions" in s:
        salida.append(l)
        salida.append(inc_lig)
        puesto_lig = True
        continue
    salida.append(l)
    # Y una molecula de ligando en [ molecules ], que es lo unico que hay que
    # tocar ahi. Hay que saber que estamos DENTRO de esa seccion: la cadena
    # "Protein_chain_A" aparece antes, en [ moleculetype ], y meter la linea
    # ahi declara un moleculetype llamado LIG y GROMACS aborta con
    # "moleculetype LIG is redefined".
    #
    # El orden importa: en [ molecules ] el ligando tiene que ir justo detras
    # de la proteina y antes del disolvente, porque esa lista tiene que seguir
    # el orden de los atomos del .gro, y el complejo es proteina + ligando.
    if s.startswith("[") and s.endswith("]"):
        en_molecules = s.strip("[] ").strip() == "molecules"
        continue
    if en_molecules and not puesto_mol and not s.startswith(";") and s:
        partes = s.split()
        if len(partes) == 2 and partes[1].isdigit():
            salida.append(f"{'LIG':<20s}1")
            puesto_mol = True

if not (puesto_tipos and puesto_lig and puesto_mol):
    sys.exit(f"no he sabido editar {top}: tipos={puesto_tipos} "
             f"ligando={puesto_lig} molecules={puesto_mol}")

top.write_text("\n".join(salida) + "\n")
print(f"topologia: {top} actualizada")
print(f"  {inc_tipos}")
print(f"  {inc_lig}")
print("  LIG 1 en [ molecules ]")
