# Superpone las cuatro conformaciones por las que pasa el ligando 453 -- la de
# la biblioteca de Pharmit, la pose de GNINA, el inicio de la dinamica y los
# 50 ns -- y renderiza las comparaciones por parejas.
# Salida: run/render/conformeros/par_*.png
#
# Documenta algo que no se ve en ninguna otra figura del repositorio: que el
# docking retuerce el enlazador del benzotriazol para que la pose encaje. En la
# biblioteca ese diedro esta en anti (-179,9); GNINA lo deja en 115,2, a 5
# grados de la cima de su barrera torsional, que el campo de fuerza cifra en
# 3,4 kcal/mol. La funcion de puntuacion no cobra eso: Vina penaliza el NUMERO
# de enlaces rotables, no el valor que toman.
#
# En ChimeraX, desde la raiz del repositorio:
#   open data/md/ligand_453/frame0.pdb
#   open data/md/ligand_453/traj_ligera.xtc structureModel #1
#   runscript scripts/md_analysis/conformeros.py
import os
import re
import sys
from collections import defaultdict, deque

import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.core.commands import run

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import exige_trayectoria  # noqa: E402

BASE = "run/render/conformeros/"       # en .gitignore
SDF = "results/pharmit/top1000_rmsd.sdf"
MOL2 = "results/docking/candidatos/docked_ligando_453.mol2"
COMPUESTO = "000-056-093"
COV = {"H": 0.38, "C": 0.77, "N": 0.75, "O": 0.73, "S": 1.02}
COLOR = {"pharmit": "#2563eb", "gnina": "#ea580c",
         "md_t0": "#0f766e", "md_t50": "#7c3aed"}
PAREJAS = [(("pharmit", "gnina"), "1_docking"),
           (("gnina", "md_t0"), "2_minimizacion"),
           (("md_t0", "md_t50"), "3_dinamica"),
           (("pharmit", "md_t50"), "4_vuelta")]

os.makedirs(BASE, exist_ok=True)


# --- Lectura de los formatos de cada etapa -----------------------------------
def lee_sdf(path, marca):
    """Extrae del SDF multimolecula el registro cuyo nombre contiene 'marca'."""
    bloque = [b for b in open(path).read().split("$$$$\n")
              if marca in b.split("\n")[0]][0]
    L = bloque.splitlines()
    na, nb = int(L[3][0:3]), int(L[3][3:6])
    els = [ln[31:34].strip() for ln in L[4:4 + na]]
    xyz = [[float(ln[0:10]), float(ln[10:20]), float(ln[20:30])]
           for ln in L[4:4 + na]]
    bonds = [(int(ln[0:3]) - 1, int(ln[3:6]) - 1) for ln in L[4 + na:4 + na + nb]]
    return els, np.array(xyz), bonds


def lee_mol2(path):
    """El mol2 de GNINA: RDKit se atraganta con sus tipos Tripos, asi que se
    leen las secciones a mano."""
    t = open(path).read()
    at = t.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>BOND")[0].strip().splitlines()
    bo = t.split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0].strip().splitlines()
    els, xyz = [], []
    for ln in at:
        f = ln.split()
        els.append(re.sub(r"[^A-Za-z]", "", f[5].split(".")[0]))
        xyz.append([float(x) for x in f[2:5]])
    bonds = [(int(l.split()[1]) - 1, int(l.split()[2]) - 1) for l in bo]
    return els, np.array(xyz), bonds


def enlaces_por_distancia(els, X):
    """Para las conformaciones de la dinamica, que vienen sin tabla de enlaces."""
    bonds = []
    for i in range(len(X)):
        d = np.linalg.norm(X - X[i], axis=1)
        for j in np.where((d > 0.1) & (d < 2.0))[0]:
            if j > i and d[j] < (COV.get(els[i], 0.8) + COV.get(els[j], 0.8)) * 1.25:
                bonds.append((i, int(j)))
    return bonds


# --- Marco rigido comun -------------------------------------------------------
def anclas(els, bonds):
    """Devuelve los indices del marco sobre el que superponer: los 6 carbonos
    del anillo central, los 2 oxigenos de metoxilo y los 2 nitrogenos amidicos.

    Las amidas NO sirven de anclaje: GNINA rota los enlaces arilo-N, asi que su
    posicion relativa cambia entre etapas. Solo el anillo y lo que cuelga
    directamente de el es rigido en las cuatro.
    """
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    pes = lambda i: [j for j in adj[i] if els[j] != "H"]   # noqa: E731

    def ciclos(n):
        out, vistos = [], set()
        for ini in range(len(els)):
            pila = [(ini, [ini])]
            while pila:
                x, cam = pila.pop()
                for y in pes(x):
                    if y == ini and len(cam) == n:
                        k = frozenset(cam)
                        if k not in vistos:
                            vistos.add(k)
                            out.append(list(cam))
                    elif y not in cam and len(cam) < n:
                        pila.append((y, cam + [y]))
        return out

    anillos6, todos = ciclos(6), ciclos(6) + ciclos(5)
    amidas = {}
    for n in range(len(els)):
        if els[n] != "N":
            continue
        for c in pes(n):
            if els[c] != "C":
                continue
            o = [x for x in pes(c) if els[x] == "O" and len(pes(x)) == 1]
            if not o:
                continue
            acilo = [x for x in pes(c) if x not in (n, o[0])][0]
            # el acilo del tiofeno esta en un anillo; el del benzotriazol es un CH2
            grupo = "tiofeno" if any(acilo in a for a in todos) else "benzotriazol"
            amidas[grupo] = (n, c, o[0])

    Nt, Nb = amidas["tiofeno"][0], amidas["benzotriazol"][0]
    ct = [x for x in pes(Nt) if any(x in a for a in anillos6)][0]
    cb = [x for x in pes(Nb) if any(x in a for a in anillos6)][0]
    central = [a for a in anillos6
               if ct in a and cb in a and all(els[i] == "C" for i in a)][0]

    dist, cola = {Nt: 0}, deque([Nt])
    while cola:
        x = cola.popleft()
        for y in pes(x):
            if y not in dist:
                dist[y] = dist[x] + 1
                cola.append(y)

    # El anillo es para-disustituido: dar la vuelta por un lado o por otro
    # cuesta lo mismo, asi que la distancia no desempata y el recorrido puede
    # salir invertido entre ficheros. Desempata la quimica: de los dos vecinos
    # de ct dentro del anillo, solo uno lleva metoxilo.
    vecinos_anillo = [y for y in pes(ct) if y in central]
    con_O = [y for y in vecinos_anillo if any(els[z] == "O" for z in pes(y))]
    inicio = (con_O[0] if len(con_O) == 1
              else min(vecinos_anillo, key=lambda y: dist.get(y, 99)))
    orden = [ct, inicio]
    while len(orden) < 6:
        orden.append([y for y in pes(orden[-1])
                      if y in central and y != orden[-2]][0])

    metox = sorted((dist.get(o, 99), o) for o in range(len(els))
                   if els[o] == "O" and len(pes(o)) == 2
                   and any(c in central for c in pes(o)))
    return orden + [o for _, o in metox[:2]] + [Nt, Nb]


def kabsch(P, Q):
    cp, cq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cq - R @ cp


# --- Las cuatro etapas --------------------------------------------------------
s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
exige_trayectoria(s)
lig = s.residues.filter(s.residues.names == "LIG")[0]
ats = list(lig.atoms)
els_md = [a.element.name for a in ats]
idx_md = np.array([a.coord_index for a in ats])
bonds_md = enlaces_por_distancia(
    els_md, s.coordset(list(s.coordset_ids)[0]).xyzs[idx_md])

etapas = {"pharmit": lee_sdf(SDF, COMPUESTO), "gnina": lee_mol2(MOL2)}
for cid, nom in ((min(s.coordset_ids), "md_t0"), (max(s.coordset_ids), "md_t50")):
    etapas[nom] = (els_md, s.coordset(cid).xyzs[idx_md].copy(), bonds_md)

ref, vista = None, None
for nom in ("pharmit", "gnina", "md_t0", "md_t50"):
    els, X, bonds = etapas[nom]
    A = anclas(els, bonds)
    if ref is None:
        ref, Y = X[A].copy(), X
        # Orientacion comun: el anillo de frente y el eje N->N horizontal, para
        # que la vista por defecto ya ensene lo que interesa.
        anillo = X[A[:6]]
        e1 = X[A[9]] - X[A[8]]
        e1 /= np.linalg.norm(e1)
        e3 = np.linalg.svd(anillo - anillo.mean(0))[2][2]
        e2 = np.cross(e3, e1)
        e2 /= np.linalg.norm(e2)
        vista = (np.vstack([np.cross(e2, e3), e2, e3]), anillo.mean(0))
    else:
        R, t = kabsch(X[A], ref)
        Y = (R @ X.T).T + t
        rms = np.sqrt((((R @ X[A].T).T + t - ref) ** 2).sum(1).mean())
        session.logger.info(f"{nom:8s} RMSD del nucleo tras superponer: {rms:.3f} A")
    Rv, c0 = vista
    Y = (Rv @ (Y - c0).T).T
    with open(f"{BASE}sup_{nom}.pdb", "w") as f:
        for i, (e, p) in enumerate(zip(els, Y)):
            f.write(f"HETATM{i + 1:5d} {e + str(i + 1):<4s} LIG A   1    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {e:>2s}\n")
        for a, b in bonds:
            f.write(f"CONECT{a + 1:5d}{b + 1:5d}\n")
        f.write("END\n")

# --- Render -------------------------------------------------------------------
run(session, "close session")
orden = ["pharmit", "gnina", "md_t0", "md_t50"]
for i, nom in enumerate(orden, start=1):
    run(session, f"open {BASE}sup_{nom}.pdb")
    run(session, f"style #{i} stick")
    run(session, f"size #{i} stickRadius 0.16")
    run(session, f"color #{i} {COLOR[nom]}")
    run(session, f"hide #{i} & H atoms")
run(session, "set bgColor white")
run(session, "material dull")
run(session, "lighting full intensity 0.9 shadows false ambientIntensity 1.05")
run(session, "graphics silhouettes true width 1.1 color #1b1f24")
run(session, "view")
run(session, "zoom 0.85")
for (a, b), nombre in PAREJAS:
    run(session, "hide #1-4 models")
    run(session, f"show #{orden.index(a) + 1},{orden.index(b) + 1} models")
    run(session, f"save {BASE}par_{nombre}.png width 1400 height 1000 supersample 3")
session.logger.info(f"Renderizadas {len(PAREJAS)} parejas en {BASE}")
