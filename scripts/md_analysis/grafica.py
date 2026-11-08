import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Por defecto, la raiz del repositorio.
d = sys.argv[1] if len(sys.argv) > 1 else "."
TABLAS = d + "/results/md/tables"
FIGURAS = d + "/results/md/figures"
fr = list(csv.DictReader(open(TABLAS + "/por_fotograma.csv")))
an = list(csv.reader(open(TABLAS + "/anillos.csv")))[1:]
t = np.array([float(x["tiempo_ns"]) for x in fr])
rmsd = np.array([float(x["rmsd_lig_A"]) for x in fr])
benzotriazol = np.array([float(x[4]) for x in an])
dimetoxifenilo = np.array([float(x[3]) for x in an])
gly = np.array([float(x[5]) for x in an])


def smooth(y, w=10):
    p = np.pad(y, (w // 2, w - 1 - w // 2), mode="edge")
    return np.convolve(p, np.ones(w) / w, mode="valid")


plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True, gridspec_kw={"height_ratios": [1, 1.2]})
for ax in (a1, a2):
    ax.axvspan(27.5, 30, color="#f2c14e", alpha=0.25, lw=0)
a1.plot(t, rmsd, color="#9aa5b1", lw=0.6)
a1.plot(t, smooth(rmsd), color="#1f2933", lw=2)
a1.set_ylabel("RMSD ligando (Å)")
a1.set_title("Ligando 453 en PARP1 — cambio de modo de unión a ~27,5 ns", loc="left", fontweight="bold")
a1.text(28.75, 6.2, "transición", ha="center", color="#b07d00")
a2.plot(t, smooth(benzotriazol), color="#e0457b", lw=2.2, label="benzotriazol")
a2.plot(t, smooth(dimetoxifenilo), color="#2aa198", lw=2.2, label="dimetoxifenilo")
a2.plot(t, smooth(gly), color="#7b8794", lw=1.5, ls="--", label="Gly863 N → aceptor más cercano")
a2.axhline(4.5, color="#cbd2d9", lw=1, ls=":")
a2.text(0.5, 4.65, "dentro del sándwich Tyr896/Tyr907", color="#7b8794", fontsize=10)
a2.set_ylabel("Distancia al centro\nTyr896/Tyr907 (Å)")
a2.set_xlabel("Tiempo (ns)")
a2.set_ylim(1, 10)
a2.legend(frameon=False, loc="upper left", ncol=3, fontsize=10)
a2.set_xlim(0, 50)
fig.tight_layout()
fig.savefig(FIGURAS + "/transicion.png", dpi=150)
