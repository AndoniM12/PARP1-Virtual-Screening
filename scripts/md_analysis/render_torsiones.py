# Renderiza el ligando solo, a 60 fps, con los tres enlaces del enlazador del
# benzotriazol resaltados: son los unicos que cambian de rotamero.
# Salida: run/render/torsiones/mol/NNNN.png
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/frame0.pdb
#   open data/md/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/render_torsiones.py
#
# NO ejecutar suaviza.py antes: la media movil amortigua las oscilaciones de
# los diedros que este video enseña.
import csv
import math
import os
import sys

import numpy as np
from chimerax.atomic import Atom, all_atomic_structures
from chimerax.core.commands import run
from chimerax.geometry import align_points

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import encuentra_amidas, vecinos  # noqa: E402

BASE = "run/render/torsiones/"         # en .gitignore
PASO = float(os.environ.get("TORS_PASO", "0.5"))   # coordsets por fotograma
ANCHO, ALTO = 1200, 900
TEMPORAL = 9999
OSCILA = 15.0                          # grados de vaiven de la camara

GRIS = (138, 150, 163, 255)            # el nucleo, que no se mueve
BRAZO = (196, 181, 253, 255)           # el benzotriazol, que si
# Un color por enlace del enlazador, en el orden en que salen del carbonilo.
COLORES = [(245, 158, 11, 255), (109, 40, 217, 255), (8, 145, 178, 255)]

os.makedirs(BASE + "mol", exist_ok=True)

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
lig = s.residues.filter(s.residues.names == "LIG")[0]
amidas = encuentra_amidas(s, lig)

# --- Nucleo, brazo y los tres enlaces ----------------------------------------
anillos = [list(r.ordered_atoms)
           for r in s.rings(cross_residue=False, all_size_threshold=6)
           if list(r.ordered_atoms)[0].residue is lig]
fusion = {a for ats in anillos if sum(x.element.name == "N" for x in ats) >= 2
          for a in ats}
for ats in anillos:
    if set(ats) & fusion and len(ats) == 6:
        fusion |= set(ats)

bt = [a for a in amidas if a[0] == "benzotriazol"][0]
cadena = [bt[2], bt[1]]                # carbonilo, C_acilo
while True:
    siguientes = [x for x in vecinos(cadena[-1]) if x is not cadena[-2]]
    if len(siguientes) != 1:
        break
    cadena.append(siguientes[0])
enlaces = [(cadena[i], cadena[i + 1]) for i in range(len(cadena) - 1)]
brazo = fusion | set(cadena[1:-1])
session.logger.info(f"enlaces del enlazador: {len(enlaces)}")

# --- Superposicion sobre el NUCLEO -------------------------------------------
# Sobre el ligando entero, el ajuste reparte el movimiento del brazo por toda
# la molecula y parece que se mueve todo. Ajustando solo el nucleo, el nucleo
# queda quieto y el brazo se ve salir, que es lo que el video cuenta.
pesados = [a for a in lig.atoms if a.element.name != "H"]
idx = np.array([a.coord_index for a in pesados if a not in brazo])
todos = np.array([a.coord_index for a in pesados])
ids = [c for c in s.coordset_ids if c != TEMPORAL]
ref = s.coordset(ids[0]).xyzs[idx]
xyz = {}
for cid in ids:
    X = s.coordset(cid).xyzs.copy()
    xyz[cid] = align_points(X[idx], ref)[0].transform_points(X)
    s.add_coordset(cid, xyz[cid])

# --- Escena -------------------------------------------------------------------
run(session, "2dlabels delete")        # el texto lo pone componer_torsiones.py
run(session, "set bgColor white")
run(session, "hide #1 cartoons")
run(session, "hide #1 atoms")
run(session, "style #1 stick")

for a in pesados:
    a.display = True
    a.draw_mode = Atom.STICK_STYLE
    a.color = BRAZO if a in brazo else GRIS

for b in lig.atoms.intra_bonds:
    b.radius = 0.26
# Los tres enlaces del enlazador, cada uno de su color y mas gruesos. Hay que
# apagar halfbond: por defecto el enlace se pinta con el color de sus atomos.
destacados = {frozenset(e) for e in enlaces}
for b in lig.atoms.intra_bonds:
    if frozenset(b.atoms) in destacados:
        i = [frozenset(e) for e in enlaces].index(frozenset(b.atoms))
        b.halfbond = False
        b.color = COLORES[i]
        b.radius = 0.46

run(session, "material dull")
run(session, "lighting full intensity 0.9 shadows true qualityOfShadows normal "
             "multiShadow 32 ambientIntensity 0.8")
run(session, "graphics silhouettes true width 1.6 color #1b1f24")
# Encuadre: con el nucleo fijo, el brazo barre varios angstrom y la vista del
# fotograma 0 se queda corta. Se encuadra sobre el fotograma mas extendido.
centro = ref.mean(axis=0)
extension = {cid: float(np.linalg.norm(xyz[cid][todos] - centro, axis=1).max())
             for cid in ids}
s.active_coordset_id = max(extension, key=extension.get)
run(session, "view #1 & :LIG")
run(session, "zoom 0.92")

# --- Fotogramas ---------------------------------------------------------------
# Los intermedios se interpolan linealmente entre coordsets consecutivos: a
# 60 fps, reproducir los 1001 originales sin mas daria un video el doble de
# rapido, no uno mas fluido.
ultimo = max(ids)
plan, f = [], 1.0
while f <= ultimo:
    plan.append(f)
    f += PASO

angulo_prev = 0.0
with open(BASE + "plan.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["salida", "coordset", "tiempo_ns"])
    for i, f in enumerate(plan):
        a = int(math.floor(f))
        b = min(a + 1, ultimo)
        x = f - a
        s.add_coordset(TEMPORAL, (1 - x) * xyz[a] + x * xyz[b])
        s.active_coordset_id = TEMPORAL
        angulo = OSCILA * math.sin(2 * math.pi * i / 900)
        run(session, f"turn y {angulo - angulo_prev:.4f} center #1:LIG", log=False)
        angulo_prev = angulo
        run(session, f"save {BASE}mol/{i:04d}.png width {ANCHO} height {ALTO} "
                     f"supersample 2", log=False)
        w.writerow([i, f"{f:.2f}", f"{(f - 1) * 0.05:.3f}"])
        out.flush()
s.active_coordset_id = ids[0]
session.logger.info(f"Renderizados {len(plan)} fotogramas en {BASE}mol/")
