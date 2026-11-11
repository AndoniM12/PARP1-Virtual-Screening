# Panel que acompaña al video: un dial por amida con el diedro omega en vivo,
# un tercer dial con la torsion CH2-CH2 del enlazador del benzotriazol, y
# debajo esa misma torsion avanzando en el tiempo.
#
# Las dos agujas de las amidas no se mueven en 50 ns; la del enlazador salta
# dos veces, en 27,75 y 30,65 ns. Esa es toda la historia del video.
# Salida: run/render/amidas/inset/NNNN.png (fondo transparente)
#
#   python scripts/md_analysis/inset_amidas.py
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
RENDER = D + "/run/render/amidas"          # en .gitignore
os.makedirs(RENDER + "/inset", exist_ok=True)

COLOR = {"tiofeno": "#c2410c", "benzotriazol": "#0f766e", "enlazador": "#6d28d9"}
PLANO = 20.0            # margen de "amida plana" alrededor de 0 y de 180
ESCALONADO = 30.0       # margen de rotamero alrededor de 60, 180 y 300
TORSION = "enlazador_CH2-CH2"


def carga(p):
    f = list(csv.DictReader(open(p)))
    return f, np.array([float(x["tiempo_ns"]) for x in f])


fil, t = carga(TABLAS + "/amidas.csv")
tor, _ = carga(TABLAS + "/torsiones.csv")
grupos = [c[len("omega_"):] for c in fil[0] if c.startswith("omega_")]
W = {g: np.array([float(x["omega_" + g]) for x in fil]) for g in grupos}
# La torsion se dibuja en [0, 360): asi los tres rotameros caen en 60, 180 y
# 300, y la linea no salta de +180 a -180 cada dos fotogramas.
T = np.array([float(x[TORSION]) for x in tor]) % 360
plan = list(csv.DictReader(open(RENDER + "/plan.csv")))

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(6.4, 3.9), dpi=120)
fig.patch.set_alpha(0)
gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1], hspace=0.55, wspace=0.25)

agujas, textos = {}, {}


def dial(col, clave, titulo, marcas, margen, color):
    ax = fig.add_subplot(gs[0, col], projection="polar")
    ax.set_facecolor((1, 1, 1, 0))
    ax.set_theta_zero_location("E")
    ax.set_ylim(0, 1.0)
    ax.set_yticks([])
    ax.set_xticks(np.radians([m for m, _ in marcas]))
    ax.set_xticklabels([e for _, e in marcas], fontsize=7, color="#52606d")
    ax.tick_params(pad=-3)
    ax.grid(color="#dde2e8", lw=0.6)
    ax.spines["polar"].set_color("#cbd2d9")
    for centro, _ in marcas:
        ax.fill_between(np.radians(np.linspace(centro - margen, centro + margen, 40)),
                        0, 1.0, color="#86c08a", alpha=0.26, lw=0)
    agujas[clave], = ax.plot([0, 0], [0, 0.92], color=color, lw=3.0,
                             solid_capstyle="round", zorder=5)
    ax.plot([0], [0], "o", color=color, ms=5, zorder=6)
    ax.set_title(titulo, color=color, fontsize=10.5, fontweight="bold", pad=13)
    textos[clave] = ax.text(0.5, -0.32, "", transform=ax.transAxes, ha="center",
                            fontsize=9.5, color="#1f2933")


MARCAS_AMIDA = [(0, "0°\ncis"), (90, "90°"), (180, "180°\ntrans"), (270, "−90°")]
MARCAS_TOR = [(60, "g+"), (180, "anti"), (300, "g−")]
for k, g in enumerate(grupos):
    dial(k, g, g, MARCAS_AMIDA, PLANO, COLOR[g])
dial(2, "enlazador", "enlazador CH₂–CH₂", MARCAS_TOR, ESCALONADO, COLOR["enlazador"])

# --- La torsion del enlazador en el tiempo ------------------------------------
ax = fig.add_subplot(gs[1, :])
ax.set_facecolor((1, 1, 1, 0))
ax.axvspan(27.5, 30, color="#f2c14e", alpha=0.35, lw=0)
for centro, _ in MARCAS_TOR:
    ax.axhspan(centro - ESCALONADO, centro + ESCALONADO, color="#86c08a",
               alpha=0.16, lw=0)



def con_cortes(x, y, umbral=120):
    """Parte la linea en los saltos grandes. Entre dos fotogramas separados
    50 ps no sabemos por donde giro la torsion: unir un g+ con un g- dibujaria
    una recta que cruza 'anti' y sugeriria un camino que no hemos medido."""
    xs, ys, pos = [], [], []
    for i in range(len(x)):
        if i and abs(y[i] - y[i - 1]) > umbral:
            xs.append(np.nan)
            ys.append(np.nan)
        pos.append(len(xs))
        xs.append(x[i])
        ys.append(y[i])
    return np.array(xs), np.array(ys), np.array(pos)


Tx, Ty, Tpos = con_cortes(t, T)
ax.plot(Tx, Ty, color=COLOR["enlazador"], lw=0.6, alpha=0.20)
pasado, = ax.plot([], [], color=COLOR["enlazador"], lw=1.7)
punto, = ax.plot([], [], "o", color=COLOR["enlazador"], ms=5.5)
ax.set_xlim(0, 50)
ax.set_ylim(0, 360)
ax.set_yticks([60, 180, 300])
ax.set_yticklabels(["g+", "anti", "g−"])
ax.set_xticks([0, 10, 20, 30, 40, 50])
ax.set_xlabel("tiempo (ns)", color="#52606d", labelpad=1)
for lado in ("top", "right"):
    ax.spines[lado].set_visible(False)
for lado in ("left", "bottom"):
    ax.spines[lado].set_color("#9aa5b1")
ax.tick_params(colors="#52606d", labelsize=8.5)
fig.subplots_adjust(left=0.075, right=0.985, top=0.87, bottom=0.115)

for p in plan:
    i = int(p["coordset"]) - 1
    for g in grupos:
        agujas[g].set_data([np.radians(W[g][i])] * 2, [0, 0.92])
        estado = "cis" if abs(W[g][i]) < 90 else "trans"
        textos[g].set_text(f"ω = {W[g][i]:+.0f}°  ({estado})")
    agujas["enlazador"].set_data([np.radians(T[i])] * 2, [0, 0.92])
    rot = "anti" if abs(T[i] - 180) < 120 else ("g+" if T[i] < 180 else "g−")
    textos["enlazador"].set_text(f"τ = {T[i]:.0f}°  ({rot})")
    pasado.set_data(Tx[:Tpos[i] + 1], Ty[:Tpos[i] + 1])
    punto.set_data([t[i]], [T[i]])
    fig.savefig(f"{RENDER}/inset/{int(p['salida']):04d}.png", transparent=True)
print(f"Generados {len(plan)} paneles en {RENDER}/inset/")
