# Puede el ligando 453 alcanzar el puente de hidrogeno con GLY863 sin tensionar
# el enlazador del benzotriazol?
#
# GLY863 es el puente de hidrogeno mas conservado del farmacoforo de consenso
# (67,74 % de los complejos de referencia, results/docking/interacciones.csv) y
# lo acepta un nitrogeno del triazol: el extremo del brazo, al otro lado del
# enlazador. Asi que la pregunta se contesta barriendo los tres diedros de ese
# enlazador con todo lo demas clavado, y mirando si alguna combinacion llega sin
# chocar.
#
# Se barre en dos geometrias:
#   - la pose de GNINA, donde el puente existe (3,28 A)
#   - el fotograma 0 de la dinamica, ya minimizado y equilibrado (5,72 A)
#
# Resultado: en la geometria de la pose de GNINA el puente es alcanzable, pero
# la combinacion sin choque mas barata que lo hace cuesta 4,23 kcal/mol de
# torsion, de un maximo posible de 7,33; relajando los tres diedros (<= 1,58) el
# triazol se queda a 9,7-9,8 A. En el fotograma 0 ya no se alcanza a ningun valor
# de torsion (minimo 4,69 A). Es decir: el puente que justificaba la pose y un
# enlazador relajado son incompatibles.
#
# De paso: GNINA eligio 113/115/32, que cuesta 6,26. Hay una via a 4,23
# (-150/-130/-130) que llega a la misma distancia, asi que la pose no solo esta
# tensionada: esta mas tensionada de lo necesario para la interaccion que hace.
#
# Limites: solo terminos de torsion propia (sin no enlazantes, sin proteina, sin
# disolvente); nucleo y proteina rigidos; rejilla de PASO grados; criterio de
# choque de esfera dura, calibrado con las propias geometrias reales, que tienen
# contactos pesado-pesado de hasta 2,92 A (un corte de 3,1 A las rechazaria).
#
# Desde la raiz del repositorio:
#   python3 scripts/md_analysis/alcance.py [paso_en_grados]
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tension_lib import (Ligando, curvas_torsionales, diedro,  # noqa: E402
                         distal, distancias_topologicas, grafo,
                         lee_mol2, lee_pdb_ligando, lee_proteina,
                         lee_topologia, rota, solo_pesados, tension)

PASO = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
ITP = "data/md/ligand_453/gromacs/ligand_453_LIG.itp"
FRAME0 = "data/md/frame0.pdb"
MOL2 = "results/docking/candidatos/docked_ligando_453.mol2"
COMPLEJO = "results/docking/candidatos_complejo/docked_ligando_453_complex.pdb"
DESPL = 661          # el complejo esta renumerado: GLY863 -> GLY202
HB = 3.5             # N...N que cuenta como puente de hidrogeno
UM_POLAR, UM_RESTO, UM_INTRA = 2.5, 2.8, 2.8

els_t, bonds_t, tor = lee_topologia(ITP)
els_h, adj_h, pesados = solo_pesados(els_t, bonds_t)
lig = Ligando(els_h, adj_h)
_, _, X_f0_todo, _ = lee_pdb_ligando(FRAME0)
rej, curvas = curvas_torsionales(lig, tor, X_f0_todo, pesados)

ROT = lig.enlazador
QUAD = [lig.quad_de(e) for e in ROT]
MUEVE = [distal(adj_h, p, q) for p, q in ROT]
ETQ = ["CO-CH2 (c_acilo)", "CH2-CH2", "CH2-triazol"]
NOMBRE_ENL = [n for n, q in lig.rotables
              if {q[1], q[2]} in [set(e) for e in ROT]]
# El nitrogeno sustituido del triazol no acepta puentes: su par libre esta en el
# sistema aromatico. Se excluye.
ACC = [i for i in lig.triazol if els_h[i] == "N" and i != lig.cadena[-1]]

print("cadena del enlazador: "
      + " - ".join(f"{els_h[i]}{i + 1}" for i in lig.cadena))
for k, (p, q) in enumerate(ROT):
    e, n = curvas[lig.rotables[[i for i, (nm, _) in enumerate(lig.rotables)
                                if {lig.rotables[i][1][1],
                                    lig.rotables[i][1][2]} == {p, q}][0]][0]]
    print(f"  {els_h[p]}{p + 1}-{els_h[q]}{q + 1:<3d} {ETQ[k]:18s} "
          f"barrera {e.max():5.2f} kcal/mol, {n} terminos, "
          f"mueve {len(MUEVE[k])} atomos pesados")
print("aceptores del triazol: " + ", ".join(f"N{i + 1}" for i in ACC))
BARRERA_TOTAL = sum(curvas[n][0].max() for n in NOMBRE_ENL)

dt = distancias_topologicas(adj_h)
ats = sorted(adj_h)
PARES = [(i, j) for a, i in enumerate(ats) for j in ats[a + 1:]
         if dt[i].get(j, 9) > 3]
PI = np.array([p[0] for p in PARES])
PJ = np.array([p[1] for p in PARES])


def energia(phis):
    return sum(tension(curvas, rej, NOMBRE_ENL[k], phis[k]) for k in range(3))


def barre(titulo, Xh, prot, gly_n):
    Pxyz = np.array([p[4] for p in prot])
    Ppol = np.array([p[3] in ("N", "O") for p in prot])
    movs = sorted(set(MUEVE[0]))
    mov_pol = np.array([els_h[i] in ("N", "O") for i in movs])
    centro = Xh[ROT[0][0]]
    cerca = np.linalg.norm(Pxyz - centro, axis=1) < (
        np.linalg.norm(Xh[movs] - centro, axis=1).max() + 6.0)
    Pc = Pxyz[cerca]
    UM = np.where(mov_pol[:, None] & Ppol[cerca][None, :], UM_POLAR, UM_RESTO)

    phis0 = [diedro(Xh, QUAD[k]) for k in range(3)]
    d0 = min(np.linalg.norm(Xh[i] - gly_n) for i in ACC)
    print(f"\n=== {titulo} ===")
    print(f"  GLY863 N -> triazol: {d0:.2f} A")
    print("  de partida: " + "  ".join(
        f"{ETQ[k].split()[0]} = {phis0[k]:7.1f}°" for k in range(3))
        + f"   E_tors = {energia(phis0):.2f} de {BARRERA_TOTAL:.2f} kcal/mol")

    grid = np.arange(-180, 180, PASO)
    combos = []
    for a2 in grid:
        X2 = rota(Xh, *ROT[2], MUEVE[2], a2 - diedro(Xh, QUAD[2]))
        for a1 in grid:
            X1 = rota(X2, *ROT[1], MUEVE[1], a1 - diedro(X2, QUAD[1]))
            for a0 in grid:
                Y = rota(X1, *ROT[0], MUEVE[0], a0 - diedro(X1, QUAD[0]))
                d = min(np.linalg.norm(Y[i] - gly_n) for i in ACC)
                combos.append((d, energia([a0, a1, a2]), (a0, a1, a2), Y))
    D = np.array([c[0] for c in combos])
    E = np.array([c[1] for c in combos])
    print(f"  {len(combos)} combinaciones (paso {PASO}°): "
          f"distancia min {D.min():.2f} A, max {D.max():.2f} A")

    def libre(k):
        d, e, ph, Y = combos[k]
        if (np.linalg.norm(Y[movs][:, None] - Pc[None], axis=2) < UM).any():
            return None
        if (np.linalg.norm(Y[PI] - Y[PJ], axis=1) < UM_INTRA).any():
            return None
        return (d, e, ph)

    cercanas = [r for r in (libre(k) for k in np.argsort(D)[:6000]) if r]
    relajadas = [r for r in (libre(k) for k in np.argsort(E)[:6000]) if r]
    if cercanas:
        cercanas.sort()
        print(f"  sin choque, las mas cercanas:")
        print(f"  {'dist':>6s}{'E_tors':>9s}" + "".join(
            f"{t.split()[0]:>12s}" for t in ETQ))
        for d, e, ph in cercanas[:6]:
            print(f"  {d:6.2f}{e:9.2f}" + "".join(f"{x:12.0f}" for x in ph))
        hb = sorted((c for c in cercanas if c[0] <= HB), key=lambda c: c[1])
        print(f"  con puente (<= {HB} A) y sin choque: {len(hb)}")
        for d, e, ph in hb[:3]:
            print(f"    las mas baratas: {e:5.2f} kcal/mol a {d:.2f} A"
                  + "  diedros " + " ".join(f"{x:5.0f}" for x in ph))
    else:
        print("  ninguna combinacion sin choque entre las mas cercanas")
    if relajadas:
        e_min = min(c[1] for c in relajadas)
        cerca_min = [c for c in relajadas if c[1] <= e_min + 0.25]
        print(f"  relajado (E_tors <= {e_min + 0.25:.2f}, {len(cerca_min)} "
              f"combinaciones sin choque): el triazol queda a "
              f"{min(c[0] for c in cerca_min):.2f}-"
              f"{max(c[0] for c in cerca_min):.2f} A de GLY863")


# --- geometria 1: la pose de GNINA en su complejo ----------------------------
els_m, X_m, bonds_m = lee_mol2(MOL2, pose=1)
e_h, a_h, pes_m = solo_pesados(els_m, bonds_m)
orden = lig.mapa_desde(e_h, a_h)
Xh_gnina = np.array([X_m[pes_m[j]] for j in orden])
prot_c = lee_proteina(COMPLEJO, excluir="UNL", desplaza_residuos=DESPL)
gly_c = [p for p in prot_c if p[0] == "GLY" and p[1] == 863 and p[2] == "N"]
if not gly_c:
    raise SystemExit(f"no encuentro GLY863 N en {COMPLEJO} con DESPL={DESPL}")
barre("pose de GNINA (donde el puente existe)", Xh_gnina, prot_c,
      np.array(gly_c[0][4]))

# --- geometria 2: el fotograma 0 de la dinamica ------------------------------
prot_f0 = lee_proteina(FRAME0, excluir="LIG")
gly_f0 = [p for p in prot_f0 if p[0] == "GLY" and p[1] == 863 and p[2] == "N"]
barre("fotograma 0 de la dinamica (ya equilibrado)",
      np.array([X_f0_todo[a] for a in pesados]), prot_f0,
      np.array(gly_f0[0][4]))
