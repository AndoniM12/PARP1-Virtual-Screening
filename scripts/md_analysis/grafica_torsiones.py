# Figura del barrido torsional: las diez torsiones del ligando en el tiempo,
# mas la actividad total, para ver de un vistazo donde estan los cambios.
# Lee results/md/tables/torsiones.csv y escribe
# results/md/figures/torsiones.png
#
#   python scripts/md_analysis/grafica_torsiones.py
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
FIGURAS = D + "/results/md/figures"
os.makedirs(FIGURAS, exist_ok=True)

VENTANA = (27.5, 30.0)          # la transicion, segun results/md/README.md
CASCADA = (27.75, 35.0)         # desde el disparo hasta que el brazo se asienta
V = 20                          # fotogramas (1 ns) a cada lado, como torsiones.py
CLAVE = "enlazador_CH2-CH2"
COLOR = {CLAVE: "#6d28d9", "tiofeno_CO-Nam": "#c2410c",
         "benzotriazol_CO-Nam": "#0f766e"}

f = list(csv.DictReader(open(TABLAS + "/torsiones.csv")))
t = np.array([float(x["tiempo_ns"]) for x in f])
nombres = [c for c in f[0] if c != "tiempo_ns"]
BRUTO = {n: np.array([float(x[n]) for x in f]) % 360 for n in nombres}


def corte(x, bins=72):
    """Angulo por el que conviene abrir el circulo para dibujarlo en plano:
    el menos visitado. Cortar siempre en 0 parte por la mitad las torsiones
    cuyo minimo esta justo ahi (los metoxilos), y la linea parece datos
    perdidos cuando en realidad no falta nada."""
    h, bordes = np.histogram(x, bins=bins, range=(0, 360))
    i = int(np.argmin(h))
    return (bordes[i] + bordes[i + 1]) / 2


CORTE = {n: corte(BRUTO[n]) for n in nombres}
X = {n: (BRUTO[n] - CORTE[n]) % 360 for n in nombres}


def media_circular(v):
    return np.degrees(np.arctan2(np.sin(np.radians(v)).mean(),
                                 np.cos(np.radians(v)).mean()))


def actividad(x):
    """Cuanto cambia la torsion en cada instante: diferencia entre la media
    circular del nanosegundo anterior y la del siguiente."""
    d = [abs((media_circular(x[i:i + V]) - media_circular(x[i - V:i]) + 180)
             % 360 - 180) for i in range(V, len(x) - V)]
    return np.array(d)


def con_cortes(x, y, umbral=120):
    """Parte la linea en los saltos. Entre dos fotogramas separados 50 ps no
    sabemos por donde giro la torsion; unir los puntos dibujaria un camino
    que no hemos medido."""
    xs, ys = [], []
    for i in range(len(x)):
        if i and abs(y[i] - y[i - 1]) > umbral:
            xs.append(np.nan)
            ys.append(np.nan)
        xs.append(x[i])
        ys.append(y[i])
    return np.array(xs), np.array(ys)


orden = sorted(nombres, key=lambda n: -actividad(X[n]).max())
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(len(orden) + 1, 1, figsize=(9.5, 10.5), sharex=True,
                         gridspec_kw={"height_ratios": [1] * len(orden) + [1.8],
                                      "hspace": 0.0})
fig.patch.set_facecolor("white")

for ax, n in zip(axes, orden):
    c = COLOR.get(n, "#9aa5b1")
    ax.axvspan(*CASCADA, color="#ede9fe", lw=0)
    ax.axvspan(*VENTANA, color="#f2c14e", alpha=0.35, lw=0)
    ax.plot(*con_cortes(t, X[n]), color=c, lw=0.9)
    ax.set_ylim(-20, 380)
    ax.set_yticks([])
    # La etiqueta va dentro del panel: con las filas pegadas, una etiqueta al
    # margen no deja claro a que traza pertenece.
    ax.text(0.006, 0.97, n, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, color=c, fontweight="bold" if n in COLOR else "normal",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.4))
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.spines["bottom"].set_color("#e3e8ee")
    ax.set_facecolor("white")
    if n == CLAVE:
        # Solo en la fila que importa se marcan los tres rotameros, porque el
        # eje vertical ya no tiene origen comun despues de cortar el circulo.
        for centro, etiqueta in ((60, "g+"), (180, "anti"), (300, "g−")):
            y = (centro - CORTE[n]) % 360
            ax.axhline(y, color=c, lw=0.5, alpha=0.3)
            ax.text(50.3, y, etiqueta, fontsize=7, color=c, va="center")

# Actividad total del ligando
ax = axes[-1]
tt = t[V:len(t) - V]
total = sum(actividad(X[n]) for n in orden)
ax.axvspan(*CASCADA, color="#ede9fe", lw=0)
ax.axvspan(*VENTANA, color="#f2c14e", alpha=0.35, lw=0)
ax.fill_between(tt, 0, total, color="#52606d", alpha=0.75, lw=0)
ax.set_ylabel("actividad\ntorsional (°)", rotation=0, ha="right", va="center",
              fontsize=8, color="#1f2933", labelpad=8)
ax.set_xlim(0, 50)
ax.set_xlabel("tiempo (ns)", color="#52606d")
ax.set_yticks([0, 150, 300])
for lado in ("top", "right"):
    ax.spines[lado].set_visible(False)
for lado in ("left", "bottom"):
    ax.spines[lado].set_color("#9aa5b1")
ax.tick_params(colors="#52606d", labelsize=8)

axes[0].set_title(
    "Las diez torsiones del ligando 453\n"
    "amarillo: la ventana de la transición · violeta: la cascada completa",
    fontsize=10, color="#1f2933", pad=10)
fig.subplots_adjust(left=0.135, right=0.955, top=0.935, bottom=0.055)
fig.savefig(FIGURAS + "/torsiones.png", dpi=160, facecolor="white")
print("escrito " + FIGURAS + "/torsiones.png")
