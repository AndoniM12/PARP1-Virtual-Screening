# Reescribe las cuatro conformaciones superpuestas que produce conformeros.py
# con nombres de atomo IDENTICOS en los cuatro modelos, y monta una escena de
# ChimeraX para girar los enlaces a mano comparandolas.
#
# Hace falta porque cada etapa viene de un formato distinto y los sup_*.pdb
# conservan el orden de su fichero de origen: el mismo atomo se llama C7 en uno
# y C19 en otro, asi que un comando `torsion` no sirve para los cuatro. Aqui se
# renombran todos al orden de la topologia de GROMACS.
#
# Salida: run/render/alcance/cuatro/{etapa}.pdb  y  girar_cuatro.cxc
#
# Desde la raiz del repositorio, con run/render/conformeros/ ya generado
# (ver results/md/ligand_453/README.md):
#   python3 scripts/md_analysis/escena_cuatro.py
#   chimerax run/render/alcance/cuatro/girar_cuatro.cxc
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tension_lib import (Ligando, diedro, grafo, lee_pdb_ligando,  # noqa: E402
                         lee_topologia, solo_pesados)

ITP = "data/md/ligand_453/gromacs/ligand_453_LIG.itp"
ORIGEN = "run/render/conformeros/"
BASE = "run/render/alcance/cuatro/"
ETAPAS = [("pharmit", "#2563eb"), ("gnina", "#ea580c"),
          ("md_t0", "#0f766e"), ("md_t50", "#7c3aed")]

if not os.path.isdir(ORIGEN):
    raise SystemExit(
        f"falta {ORIGEN}: generalo antes con conformeros.py en ChimeraX\n"
        "  open data/md/ligand_453/frame0.pdb\n"
        "  open data/md/ligand_453/traj_ligera.xtc structureModel #1\n"
        "  runscript scripts/md_analysis/conformeros.py")
os.makedirs(BASE, exist_ok=True)

els_t, bonds_t, _ = lee_topologia(ITP)
els_h, adj_h, pesados = solo_pesados(els_t, bonds_t)
lig = Ligando(els_h, adj_h)
print(f"enlaces rotables: {len(lig.rotables)}")

# --- reescritura de los cuatro PDB -------------------------------------------
coords = {}
for etapa, _ in ETAPAS:
    _, els, X, bonds = lee_pdb_ligando(f"{ORIGEN}sup_{etapa}.pdb")
    e_h, a_h, pes = solo_pesados(els, bonds)
    orden = lig.mapa_desde(e_h, a_h)          # topologia -> fichero
    Y = np.array([X[pes[j]] for j in orden])
    coords[etapa] = Y
    with open(f"{BASE}{etapa}.pdb", "w") as f:
        for i in range(len(els_h)):
            p = Y[i]
            f.write(f"HETATM{i + 1:5d} {lig.nombre[i]:<4s} LIG A   1    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00"
                    f"          {els_h[i]:>2s}\n")
        for a in sorted(adj_h):
            for b in sorted(adj_h[a]):
                if a < b:
                    f.write(f"CONECT{a + 1:5d}{b + 1:5d}\n")
        f.write("END\n")
    print(f"  escrito {BASE}{etapa}.pdb")

# --- el .cxc ------------------------------------------------------------------
ABS = os.path.abspath(BASE)
spec = {nom: "/A:1@" + ",".join(lig.nombre[i] for i in q)
        for nom, q in lig.rotables}
enl = [nom for nom, q in lig.rotables
       if {q[1], q[2]} in [set(e) for e in lig.enlazador]]
acc = [lig.nombre[i] for i in lig.triazol
       if els_h[i] == "N" and i != lig.cadena[-1]]

resets = []
for k, (etapa, _) in enumerate(ETAPAS, 1):
    resets.append(f"alias reset{k} " + " ; ".join(
        f"torsion #{k}{spec[nom]} {diedro(coords[etapa], q):.1f}"
        for nom, q in lig.rotables))
resets.append("alias reset " + " ; ".join(f"reset{k}" for k in range(1, 5)))

cxc = f"""# Las cuatro conformaciones del ligando 453, superpuestas por el anillo
# central, con los atomos nombrados igual en los cuatro modelos para poder
# girarlas y compararlas con los mismos comandos.
#
#   chimerax {ABS}/girar_cuatro.cxc
#
# El raton derecho queda en modo "bond rotation": pincha un enlace y arrastra.
# Mueve el lado pequeño del enlace, asi que el nucleo no se descoloca; con
# Mayus pulsada mueve el lado grande y SI rompe la superposicion.
{chr(10).join(f"open {ABS}/{e}.pdb" for e, _ in ETAPAS)}
{chr(10).join(f"rename #{k} {e}" for k, (e, _) in enumerate(ETAPAS, 1))}
style #1-4 stick
size #1-4 stickRadius 0.16
{chr(10).join(f"color #{k} {c}" for k, (_, c) in enumerate(ETAPAS, 1))}
set bgColor white
material dull
lighting full intensity 0.9 shadows false ambientIntensity 1.05
graphics silhouettes true width 1.1 color #1b1f24
view
zoom 0.85

# Distancia en vivo entre los aceptores del triazol de la biblioteca (#1) y los
# de los 50 ns (#4): mide cuanto falta para cuadrar un brazo con el otro.
{chr(10).join(f"distance #1/A:1@{n} #4/A:1@{n}" for n in acc)}
distance style radius 0.03 color #7c3aed dashes 6
distance style decimalPlaces 2

# Atajos. $1 es el numero de modelo: "enlazador 1", "todas 4".
alias enlazador {" ; ".join(f"torsion #$1{spec[n]}" for n in enl)}
alias todas {" ; ".join(f"torsion #$1{spec[n]}" for n, _ in lig.rotables)}
alias solo hide #1-4 models ; show #$1 models
alias cuatro show #1-4 models
{chr(10).join(resets)}

# Un reset por modelo porque el modo de raton "bond rotation" NO entra en la
# pila de undo (mouse_rot.py no registra ninguna accion de undo), asi que undo
# solo deshace los `torsion` escritos a mano. Cada arrastre si queda escrito en
# el log como un comando torsion, por si se prefiere revertirlos uno a uno.

# Los enlaces rotables y los cuatro atomos que definen cada torsion:
{chr(10).join(f"#   {nom:26s} {spec[nom]}" for nom, _ in lig.rotables)}
#
# Para fijar un valor:  torsion #1{spec['enlazador_CH2-CH2']} -171
# Para leer uno:        torsion #1{spec['enlazador_CH2-CH2']}

mousemode rightMode "bond rotation"
"""
open(BASE + "girar_cuatro.cxc", "w").write(cxc)
print(f"escrito {BASE}girar_cuatro.cxc")
print("\ndiedros de partida:")
print(f"  {'enlace':26s}" + "".join(f"{e:>10s}" for e, _ in ETAPAS))
for nom, q in lig.rotables:
    print(f"  {nom:26s}"
          + "".join(f"{diedro(coords[e], q):10.1f}" for e, _ in ETAPAS))
