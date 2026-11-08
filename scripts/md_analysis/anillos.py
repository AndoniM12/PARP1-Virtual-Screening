import csv
import numpy as np
from chimerax.atomic import all_atomic_structures

# Rutas relativas a la raiz del repositorio. ChimeraX debe
# arrancarse desde ahi, o usar el comando `cd` de ChimeraX.
OUT = "results/md/tables/"
s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
lig = s.residues.filter(s.residues.names == "LIG")[0]
res = {r.number: r for r in s.residues}


def ring_atoms(r):
    return [a for a in r.atoms if a.name in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")]


tyr = ring_atoms(res[896]) + ring_atoms(res[907])
gly_n = res[863].find_atom("N")
rings = []
for ring in s.rings(cross_residue=False, all_size_threshold=6):
    ats = [a for a in ring.ordered_atoms]
    if ats[0].residue is lig:
        el = "".join(sorted(a.element.name for a in ats))
        rings.append((f"{len(ats)}-{el}", ats))
lig_heavy = [a for a in lig.atoms if a.element.name != "H"]
acc = [a for a in lig_heavy if a.element.name in ("O", "N")]
rows = []
for cid in s.coordset_ids:
    s.active_coordset_id = cid
    sandwich = np.mean([a.scene_coord for a in tyr], axis=0)
    row = [round((cid - 1) * 0.05, 2)]
    row += [round(float(np.linalg.norm(np.mean([a.scene_coord for a in ats], axis=0) - sandwich)), 2) for _, ats in rings]
    row.append(round(float(min(np.linalg.norm(a.scene_coord - gly_n.scene_coord) for a in acc)), 2))
    rows.append(row)
s.active_coordset_id = 1
with open(OUT + "anillos.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["tiempo_ns"] + [n for n, _ in rings] + ["gly863N_aceptor_min"])
    w.writerows(rows)
session.logger.info("rings: " + ", ".join(n for n, _ in rings))
