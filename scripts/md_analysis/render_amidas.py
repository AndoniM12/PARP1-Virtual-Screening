# Renderiza el ligando solo, a lo largo de los 50 ns, con las dos amidas
# secundarias resaltadas. Salida: run/render/amidas/mol/NNNN.png
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/frame0.pdb
#   open data/md/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/render_amidas.py
#
# NO ejecutar suaviza.py antes: la media movil amortigua precisamente las
# oscilaciones del diedro que este video quiere enseñar.
import csv
import os
import sys

import numpy as np
from chimerax.atomic import Atom, all_atomic_structures
from chimerax.core.commands import run
from chimerax.geometry import align_points

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import encuentra_amidas, exige_trayectoria, vecinos  # noqa: E402

BASE = "run/render/amidas/"            # en .gitignore
PASO = int(os.environ.get("AMIDAS_PASO", "1"))   # >1 para pruebas rapidas
ANCHO, ALTO = 1200, 900
VUELTAS = 1.0                          # giro completo a lo largo del video

GRIS = (138, 150, 163, 255)
COLOR = {"tiofeno": (194, 65, 12, 255), "benzotriazol": (15, 118, 110, 255)}
VIOLETA = (109, 40, 217, 255)          # el enlazador, a juego con su dial

os.makedirs(BASE + "mol", exist_ok=True)

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
exige_trayectoria(s)   # falla antes de sobrescribir nada
lig = s.residues.filter(s.residues.names == "LIG")[0]
amidas = encuentra_amidas(s, lig)
session.logger.info(f"amidas encontradas: {len(amidas)}  ({', '.join(a[0] for a in amidas)})")

# --- Superposicion sobre el ligando -----------------------------------------
# Sin esto el ligando se pasea y da vueltas dentro de la caja de simulacion y el
# video es ilegible. Superponer cada fotograma sobre el primero quita el
# movimiento de solido rigido y deja solo el movimiento interno, que es el que
# interesa. Los diedros no se tocan: una rotacion rigida no los cambia.
ids = list(s.coordset_ids)
pesados = [a for a in lig.atoms if a.element.name != "H"]
idx = np.array([a.coord_index for a in pesados])
ref = s.coordset(ids[0]).xyzs[idx]
for cid in ids:
    xyz = s.coordset(cid).xyzs.copy()
    tf, _ = align_points(xyz[idx], ref)
    s.add_coordset(cid, tf.transform_points(xyz))

# --- Escena ------------------------------------------------------------------
run(session, "2dlabels delete")        # el texto lo pone componer_amidas.py
run(session, "set bgColor white")
run(session, "hide #1 cartoons")
run(session, "hide #1 atoms")
run(session, "style #1 stick")

# Se muestra el ligando sin hidrogenos, salvo los dos N-H de las amidas: son los
# que dejan ver de un vistazo si el grupo esta plano o retorcido.
visibles = [a for a in lig.atoms if a.element.name != "H"]
for *_, h_n in amidas:
    if h_n is not None:
        visibles.append(h_n)
for a in visibles:
    a.display = True
    a.draw_mode = Atom.STICK_STYLE
    a.color = GRIS

# Las seis posiciones de cada amida, en su color. El resto del ligando queda en
# gris: asi el ojo va solo a los dos grupos que se estan midiendo.
destacados = set()
for grupo, c_acilo, c, o, n, c_amino, h_n in amidas:
    for a in (c_acilo, c, o, n, c_amino, h_n):
        if a is None:
            continue
        a.color = COLOR[grupo]
        destacados.add(a)

# El enlazador propanamida del benzotriazol: los dos CH2 cuya torsion CH2-CH2
# cambia de rotamero en 27,75 ns, antes de que los anillos se reacomoden. Va en
# violeta, el mismo color que su dial en el panel. El primero de los dos es
# tambien el sustituyente acilo de la amida, asi que se le quita el verde: el
# verde queda para el plano de la amida y el violeta para el enlazador.
bt = [a for a in amidas if a[0] == "benzotriazol"][0]
cadena = [bt[2], bt[1]]                # carbonilo, C_acilo
while True:
    siguientes = [x for x in vecinos(cadena[-1]) if x is not cadena[-2]]
    if len(siguientes) != 1:
        break
    cadena.append(siguientes[0])
for a in cadena[1:3]:
    a.color = VIOLETA
    destacados.add(a)

# En modo stick el grosor lo manda el enlace, no el atomo: los enlaces internos
# de cada amida y del enlazador van mas gruesos, el resto del esqueleto fino.
for b in lig.atoms.intra_bonds:
    a1, a2 = b.atoms
    b.radius = 0.26 + 0.16 * (a1 in destacados and a2 in destacados)

run(session, "material dull")
run(session, "lighting full intensity 0.9 shadows true qualityOfShadows finer "
             "multiShadow 128 ambientIntensity 0.8")
run(session, "graphics silhouettes true width 1.6 color #1b1f24")
run(session, "view #1 & :LIG")
run(session, "zoom 1.30")

# --- Fotogramas ---------------------------------------------------------------
plan = ids[::PASO]
giro = 360.0 * VUELTAS / max(len(plan), 1)
with open(BASE + "plan.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["salida", "coordset", "tiempo_ns"])
    for i, cid in enumerate(plan):
        s.active_coordset_id = cid
        run(session, f"turn y {giro:.4f} center #1:LIG", log=False)
        run(session, f"save {BASE}mol/{i:04d}.png width {ANCHO} height {ALTO} "
                     f"supersample 2", log=False)
        w.writerow([i, cid, f"{(cid - 1) * 0.05:.3f}"])
        out.flush()
s.active_coordset_id = ids[0]
session.logger.info(f"Renderizados {len(plan)} fotogramas en {BASE}mol/")
