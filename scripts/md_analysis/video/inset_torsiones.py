# Panel del video de las torsiones: un dial por cada enlace del enlazador del
# benzotriazol, y debajo la consecuencia -- cuanto se ha desplazado el brazo
# respecto al nucleo.
# Salida: run/render/torsiones/inset/NNNN.png (fondo transparente)
#
#   python scripts/md_analysis/inset_torsiones.py
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
RENDER = D + "/run/render/torsiones"       # en .gitignore
os.makedirs(RENDER + "/inset", exist_ok=True)

VENTANA = (27.5, 30.0)
MARGEN = 30.0                              # margen de rotamero, en grados
# En el orden en que salen del carbonilo, con el color de su enlace en el 3D.
TORSIONES = [("benzotriazol_CO-CH2", "C(=O)–CH₂", "#f59e0b"),
             ("enlazador_CH2-CH2", "CH₂–CH₂", "#6d28d9"),
             ("enlazador_CH2-triazol", "CH₂–N", "#0891b2")]


def carga(p):
    f = list(csv.DictReader(open(p)))
    return f, np.array([float(x["tiempo_ns"]) for x in f])


tor, t = carga(TABLAS + "/torsiones.csv")
fra, _ = carga(TABLAS + "/fragmentos.csv")
X = {c: np.array([float(x[c]) for x in tor]) % 360 for c, _, _ in TORSIONES}
BRAZO = np.array([float(x["brazo_interno"]) for x in fra])
plan = list(csv.DictReader(open(RENDER + "/plan.csv")))

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(6.4, 3.9), dpi=120)
fig.patch.set_alpha(0)
gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1], hspace=0.55, wspace=0.25)

agujas, textos = {}, {}
for k, (col, etiqueta, color) in enumerate(TORSIONES):
    ax = fig.add_subplot(gs[0, k], projection="polar")
    ax.set_facecolor((1, 1, 1, 0))
    ax.set_theta_zero_location("E")
    ax.set_ylim(0, 1.0)
    ax.set_yticks([])
    ax.set_xticks(np.radians([60, 180, 300]))
    ax.set_xticklabels(["g+", "anti", "g−"], fontsize=7.5, color="#52606d")
    ax.tick_params(pad=-3)
    ax.grid(color="#dde2e8", lw=0.6)
    ax.spines["polar"].set_color("#cbd2d9")
    for centro in (60, 180, 300):
        ax.fill_between(np.radians(np.linspace(centro - MARGEN, centro + MARGEN, 40)),
                        0, 1.0, color="#86c08a", alpha=0.24, lw=0)
    agujas[col], = ax.plot([0, 0], [0, 0.92], color=color, lw=3.0,
                           solid_capstyle="round", zorder=5)
    ax.plot([0], [0], "o", color=color, ms=5, zorder=6)
    ax.set_title(etiqueta, color=color, fontsize=11, fontweight="bold", pad=13)
    textos[col] = ax.text(0.5, -0.32, "", transform=ax.transAxes, ha="center",
                          fontsize=9.5, color="#1f2933")

# --- La consecuencia: cuanto se ha ido el brazo -------------------------------
ax = fig.add_subplot(gs[1, :])
ax.set_facecolor((1, 1, 1, 0))
ax.axvspan(*VENTANA, color="#f2c14e", alpha=0.35, lw=0)
ax.plot(t, BRAZO, color="#6d28d9", lw=0.6, alpha=0.20)
pasado, = ax.plot([], [], color="#6d28d9", lw=1.8)
punto, = ax.plot([], [], "o", color="#6d28d9", ms=5.5)
ax.set_xlim(0, 50)
ax.set_ylim(0, 9)
ax.set_yticks([0, 3, 6, 9])
ax.set_xticks([0, 10, 20, 30, 40, 50])
ax.set_xlabel("tiempo (ns)", color="#52606d", labelpad=1)
ax.set_ylabel("brazo (Å)", color="#52606d", fontsize=9)
for lado in ("top", "right"):
    ax.spines[lado].set_visible(False)
for lado in ("left", "bottom"):
    ax.spines[lado].set_color("#9aa5b1")
ax.tick_params(colors="#52606d", labelsize=8.5)
fig.subplots_adjust(left=0.085, right=0.985, top=0.87, bottom=0.115)

for p in plan:
    # El coordset es fraccionario (hay fotogramas interpolados). Para las
    # agujas se toma el fotograma mas cercano: interpolar un angulo entre dos
    # valores que dan la vuelta por 360 inventaria posiciones intermedias que
    # no existen. La curva del brazo si se interpola, que no da vueltas.
    i = min(int(round(float(p["coordset"]))) - 1, len(t) - 1)
    tt = float(p["tiempo_ns"])
    for col, _, _ in TORSIONES:
        w = X[col][i]
        agujas[col].set_data([np.radians(w)] * 2, [0, 0.92])
        rot = "anti" if abs(w - 180) < 60 else ("g+" if w < 180 else "g−")
        textos[col].set_text(f"{w:.0f}°  ({rot})")
    y = float(np.interp(tt, t, BRAZO))
    j = int(np.searchsorted(t, tt, side="right"))
    pasado.set_data(np.append(t[:j], tt), np.append(BRAZO[:j], y))
    punto.set_data([tt], [y])
    fig.savefig(f"{RENDER}/inset/{int(p['salida']):04d}.png", transparent=True)
print(f"Generados {len(plan)} paneles en {RENDER}/inset/")
