# Tension torsional de los diez enlaces rotables del ligando 453, en las cuatro
# etapas por las que pasa y a lo largo de la trayectoria.
#
# Contesta cual de los enlaces carga con la tension que el docking introduce y
# cuando se suelta. Resultado: el docking anade 17,3 kcal/mol, el 92 % en tres
# enlaces, y el mayor con diferencia es benzotriazol_fenilo-Nam, que queda a
# 10,82 de una barrera de 10,87 -- practicamente en la cima. La minimizacion
# suelta 12,0 de esos 17,3; los 50 ns de dinamica no sueltan nada (de media la
# tension sube de 21,6 a 24,3), asi que el cambio de modo de union de los 28 ns
# no lo mueve la relajacion interna.
#
# El cero de cada enlace es su propio minimo torsional, asi que la cifra es
# "cuanto por encima de su optimo esta este enlace". Solo hay terminos de
# torsion propia: ni no enlazantes, ni proteina, ni disolvente. No es una
# energia de union.
#
# Salida: results/md/tables/tension_torsional.csv
#
# Desde la raiz del repositorio, con run/render/alcance/cuatro/ ya generado:
#   python3 scripts/md_analysis/escena_cuatro.py
#   python3 scripts/md_analysis/tension.py
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tension_lib import (Ligando, curvas_torsionales, diedro,  # noqa: E402
                         lee_pdb_ligando, lee_topologia, solo_pesados, tension)

ITP = "data/md/ligand_453/gromacs/ligand_453_LIG.itp"
FRAME0 = "data/md/frame0.pdb"
CUATRO = "run/render/alcance/cuatro/"
TORSIONES = "results/md/tables/torsiones.csv"
OUT = "results/md/tables/tension_torsional.csv"
ETAPAS = ["pharmit", "gnina", "md_t0", "md_t50"]
VENTANA = (27.5, 35.0)        # transicion y cola del reacomodo

if not os.path.isdir(CUATRO):
    raise SystemExit(f"falta {CUATRO}: generalo con escena_cuatro.py")

els_t, bonds_t, tor = lee_topologia(ITP)
els_h, adj_h, pesados = solo_pesados(els_t, bonds_t)
lig = Ligando(els_h, adj_h)

# Las fases de los terminos se miden sobre una geometria CON hidrogenos.
_, _, X_f0, _ = lee_pdb_ligando(FRAME0)
if len(X_f0) != len(els_t):
    raise SystemExit(f"frame0 tiene {len(X_f0)} atomos de ligando y la "
                     f"topologia {len(els_t)}: no se pueden emparejar")
rej, curvas = curvas_torsionales(lig, tor, X_f0, pesados)

# --- las cuatro etapas --------------------------------------------------------
coords = {}
for e in ETAPAS:
    nom, _, X, _ = lee_pdb_ligando(f"{CUATRO}{e}.pdb")
    idx = {n: k for k, n in enumerate(nom)}
    coords[e] = np.array([X[idx[lig.nombre[i]]] for i in range(len(els_h))])

# --- medias de la trayectoria -------------------------------------------------
# torsiones.csv elige el sustituyente de referencia por indice de atomo; aqui se
# elige por etiqueta canonica, asi que cada columna esta desplazada una
# constante. El desfase se mide comparando su fotograma 0 con la geometria
# md_t0 y se aplica a toda la columna.
filas_csv = list(csv.DictReader(open(TORSIONES)))
t = np.array([float(x["tiempo_ns"]) for x in filas_csv])
antes, despues = t < VENTANA[0], t >= VENTANA[1]

filas = []
for nombre, quad in lig.rotables:
    e_curva, n_terms = curvas[nombre]
    phis = [diedro(coords[et], quad) for et in ETAPAS]
    ten = [tension(curvas, rej, nombre, p) for p in phis]
    if nombre in filas_csv[0]:
        v = np.array([float(x[nombre]) for x in filas_csv])
        desfase = (diedro(coords["md_t0"], quad) - v[0] + 180) % 360 - 180
        serie = np.array([tension(curvas, rej, nombre, p + desfase) for p in v])
        medias = [float(serie[antes].mean()), float(serie[despues].mean())]
    else:
        print(f"aviso: {nombre} no esta en {TORSIONES}; sin medias")
        medias = [float("nan"), float("nan")]
    filas.append({
        "enlace": nombre,
        "barrera_kcal": round(float(e_curva.max()), 2),
        "n_terminos": n_terms,
        **{f"phi_{et}": round(p, 1) for et, p in zip(ETAPAS, phis)},
        **{f"tension_{et}": round(x, 2) for et, x in zip(ETAPAS, ten)},
        "media_0_27.5ns": round(medias[0], 2),
        "media_35_50ns": round(medias[1], 2),
        "pone_docking": round(ten[1] - ten[0], 2),
        "suelta_minimizacion": round(ten[1] - ten[2], 2),
        "suelta_dinamica": round(medias[0] - medias[1], 2),
    })

filas.sort(key=lambda f: -f["tension_gnina"])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(filas[0]))
    w.writeheader()
    w.writerows(filas)
print(f"escrito {OUT}")

# --- resumen ------------------------------------------------------------------
print()
print("TENSION TORSIONAL (kcal/mol sobre el minimo de cada enlace)")
print(f"{'enlace':26s}{'barrera':>8s}" + "".join(f"{e:>9s}" for e in ETAPAS)
      + f"{'pone dock':>11s}{'suelta min':>12s}")
for r in filas:
    print(f"{r['enlace']:26s}{r['barrera_kcal']:8.2f}"
          + "".join(f"{r['tension_' + e]:9.2f}" for e in ETAPAS)
          + f"{r['pone_docking']:11.2f}{r['suelta_minimizacion']:12.2f}")
tot = {e: sum(r["tension_" + e] for r in filas) for e in ETAPAS}
print(f"{'TOTAL':26s}{'':8s}" + "".join(f"{tot[e]:9.2f}" for e in ETAPAS)
      + f"{tot['gnina'] - tot['pharmit']:11.2f}"
      + f"{tot['gnina'] - tot['md_t0']:12.2f}")
print()
print("MEDIA EN LA TRAYECTORIA (kcal/mol)")
print(f"{'enlace':26s}{'0-27.5 ns':>11s}{'35-50 ns':>11s}{'suelta':>9s}")
for r in sorted(filas, key=lambda x: -x["suelta_dinamica"]):
    print(f"{r['enlace']:26s}{r['media_0_27.5ns']:11.2f}"
          f"{r['media_35_50ns']:11.2f}{r['suelta_dinamica']:9.2f}")
m1 = sum(r["media_0_27.5ns"] for r in filas)
m2 = sum(r["media_35_50ns"] for r in filas)
print(f"{'TOTAL':26s}{m1:11.2f}{m2:11.2f}{m1 - m2:9.2f}")
print()
print("Los fotogramas sueltos enganan: fenilo-Ometoxi2 da 0,13 en el fotograma")
print("de los 50 ns y 3,12 de media en ese tramo. Usa las medias donde las hay;")
print("pharmit y gnina son estructuras unicas y ahi no hay media posible.")
