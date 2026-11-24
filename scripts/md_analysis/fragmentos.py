# Separa el ligando en nucleo y brazo del benzotriazol, y mide cuanto se mueve
# cada parte: por dentro (deformacion) y dentro del bolsillo (desplazamiento).
#
# Contesta que consecuencia tiene el cambio de rotamero del enlazador que
# detecta torsiones.py. Resultado: el nucleo no cambia de forma en 50 ns
# (0,39 -> 0,42 A), todo el cambio conformacional es el brazo (1,92 -> 7,06 A),
# y el brazo arrastra al nucleo 3,3 A dentro del bolsillo.
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/ligand_453/frame0.pdb
#   open data/md/ligand_453/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/fragmentos.py
import csv
import os
import sys

import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.geometry import align_points

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import encuentra_amidas, exige_trayectoria, vecinos  # noqa: E402

# Ligando y paso entre fotogramas como argumentos:
#   runscript scripts/md_analysis/fragmentos.py 170 0.1
# Sin argumentos, el 453 y sus 50 ps. Lee lo que este cargado en #1 y escribe
# en results/md/ligand_<N>/tables/.
LIGANDO = sys.argv[1] if len(sys.argv) > 1 else "453"
DT_NS = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
OUT = f"results/md/ligand_{LIGANDO}/tables/"
os.makedirs(OUT, exist_ok=True)

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
exige_trayectoria(s)   # falla antes de sobrescribir nada
lig = s.residues.filter(s.residues.names == "LIG")[0]
amidas = encuentra_amidas(s, lig)


def brazo_benzotriazol():
    """El benzotriazol fusionado mas los CH2 que lo unen al resto."""
    anillos = [list(r.ordered_atoms)
               for r in s.rings(cross_residue=False, all_size_threshold=6)
               if list(r.ordered_atoms)[0].residue is lig]
    triazol = [a for ats in anillos
               if sum(x.element.name == "N" for x in ats) >= 2 for a in ats]
    fusion = set(triazol)
    for ats in anillos:                      # el benceno fusionado comparte
        if set(ats) & fusion and len(ats) == 6:
            fusion |= set(ats)
    bt = [a for a in amidas if a[0] == "benzotriazol"][0]
    cadena = [bt[2], bt[1]]                  # carbonilo, C_acilo
    while True:
        siguientes = [x for x in vecinos(cadena[-1]) if x is not cadena[-2]]
        if len(siguientes) != 1:
            break
        cadena.append(siguientes[0])
    return fusion | set(cadena[1:-1])


brazo = brazo_benzotriazol()
pesados = [a for a in lig.atoms if a.element.name != "H"]
iB = np.array([a.coord_index for a in pesados if a in brazo])
iN = np.array([a.coord_index for a in pesados if a not in brazo])
iT = np.array([a.coord_index for a in pesados])
session.logger.info(f"nucleo {len(iN)} atomos | brazo {len(iB)} atomos")


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


ids = list(s.coordset_ids)
X0 = s.coordset(ids[0]).xyzs
filas = []
for cid in ids:
    X = s.coordset(cid).xyzs
    # Directo, sin superponer: la trayectoria ya viene ajustada al backbone,
    # asi que esto es el movimiento dentro del bolsillo.
    pose_n, pose_b = rmsd(X[iN], X0[iN]), rmsd(X[iB], X0[iB])
    # Superponiendo el nucleo: aisla el movimiento interno del brazo respecto
    # al resto de la molecula, sin el movimiento del ligando entero.
    Y = align_points(X[iN], X0[iN])[0].transform_points(X)
    int_n, int_b = rmsd(Y[iN], X0[iN]), rmsd(Y[iB], X0[iB])
    # Superponiendo todo el ligando: cuanto cambia de forma en conjunto.
    Z = align_points(X[iT], X0[iT])[0].transform_points(X)
    filas.append([round((cid - 1) * DT_NS, 2), round(int_n, 2), round(int_b, 2),
                  round(rmsd(Z[iT], X0[iT]), 2), round(pose_n, 2), round(pose_b, 2)])
s.active_coordset_id = ids[0]

with open(OUT + "fragmentos.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["tiempo_ns", "nucleo_interno", "brazo_interno",
                 "ligando_interno", "nucleo_pose", "brazo_pose"])
    wr.writerows(filas)

A = np.array(filas)
t = A[:, 0]
session.logger.info(f"{'tramo':12s} " + " ".join(f"{c:>15s}" for c in
                    ["nucleo_interno", "brazo_interno", "nucleo_pose", "brazo_pose"]))
for nombre, m in (("0-27.5", t < 27.5), ("30-35", (t > 30) & (t <= 35)),
                  ("35-50", t > 35)):
    session.logger.info(f"{nombre:12s} " + " ".join(
        f"{A[m, k].mean():15.2f}" for k in (1, 2, 4, 5)))
