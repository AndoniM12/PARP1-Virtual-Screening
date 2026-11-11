# Genera el gráfico de RMSD en miniatura que avanza sincronizado con cada fotograma del vídeo.
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Por defecto, la raiz del repositorio.
D = sys.argv[1] if len(sys.argv) > 1 else "."
TABLAS = D + "/results/md/tables"
RENDER = D + "/run/render"          # en .gitignore
os.makedirs(RENDER + "/inset", exist_ok=True)
fr = list(csv.DictReader(open(TABLAS + "/por_fotograma.csv")))
t = np.array([float(x["tiempo_ns"]) for x in fr])
r = np.array([float(x["rmsd_lig_A"]) for x in fr])
w = 10
rs = np.convolve(np.pad(r, (w // 2, w - 1 - w // 2), mode="edge"), np.ones(w) / w, mode="valid")
plan = list(csv.DictReader(open(TABLAS + "/plan_render_60fps.csv")))

plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(4.6, 2.3), dpi=120)
fig.patch.set_alpha(0)
ax.set_facecolor((1, 1, 1, 0.92))
ax.axvspan(27.5, 30, color="#f2c14e", alpha=0.3, lw=0)
ax.plot(t, rs, color="#cbd2d9", lw=1.5)
pasado, = ax.plot([], [], color="#c2410c", lw=2.2)
punto, = ax.plot([], [], "o", color="#c2410c", ms=6)
ax.set_xlim(0, 50)
ax.set_ylim(0, 7)
ax.set_xlabel("tiempo (ns)", color="#52606d")
ax.set_ylabel("RMSD ligando (Å)", color="#52606d")
ax.set_xticks([0, 10, 20, 30, 40, 50])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(colors="#52606d")
fig.tight_layout()

for p in plan:
    tt = float(p["tiempo_ns"])
    i = int(np.searchsorted(t, tt, side="right"))
    yy = float(np.interp(tt, t, rs))
    pasado.set_data(np.append(t[:i], tt), np.append(rs[:i], yy))
    punto.set_data([tt], [yy])
    fig.savefig(f"{RENDER}/inset/{int(p['salida']):04d}.png", transparent=True)
