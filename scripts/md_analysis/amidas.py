# Mide el diedro omega (C_acilo-C-N-C_amino) de las dos amidas secundarias del
# ligando a lo largo de la trayectoria, y su desviacion de la planaridad.
#
# Contesta el punto 1 de docs/limitations.md: si alguna de las dos amidas gira
# durante la simulacion, y en particular en la ventana de la transicion.
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/frame0.pdb
#   open data/md/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/amidas.py
import csv
import os
import sys

import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.geometry import dihedral

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import desviacion_plano, encuentra_amidas  # noqa: E402

OUT = "results/md/tables/"

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
lig = s.residues.filter(s.residues.names == "LIG")[0]
amidas = encuentra_amidas(s, lig)
session.logger.info(f"amidas encontradas: {len(amidas)}")

filas = []
for cid in s.coordset_ids:
    s.active_coordset_id = cid
    fila = [round((cid - 1) * 0.05, 2)]
    for _, c_acilo, c, o, n, c_amino, _ in amidas:
        w = dihedral(c_acilo.coord, c.coord, n.coord, c_amino.coord)
        fila += [round(w, 1), round(desviacion_plano(w), 1)]
    filas.append(fila)
s.active_coordset_id = 1

cabecera = ["tiempo_ns"]
for grupo, *_ in amidas:
    cabecera += [f"omega_{grupo}", f"desv_plano_{grupo}"]
with open(OUT + "amidas.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(cabecera)
    wr.writerows(filas)

# Resumen: estado de cada amida, y que pasa en la ventana de la transicion
t = np.array([f[0] for f in filas])
vent = (t >= 27.5) & (t <= 30.0)
for k, (grupo, *_) in enumerate(amidas):
    w = np.array([f[1 + 2 * k] for f in filas])
    d = np.array([f[2 + 2 * k] for f in filas])
    cis = np.abs(w) < 90
    session.logger.info(
        f"{grupo:13s} cis {100 * cis.mean():5.1f} % | "
        f"desv_plano mediana {np.median(d):4.1f} max {d.max():5.1f} | "
        f"en 27.5-30 ns: cis {100 * cis[vent].mean():5.1f} %, "
        f"omega {w[vent].min():7.1f} a {w[vent].max():7.1f}"
    )
