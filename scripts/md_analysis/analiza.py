import csv
import os
import sys
from collections import Counter
import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.hbonds import find_hbonds

# Rutas relativas a la raiz del repositorio. ChimeraX debe
# arrancarse desde ahi, o usar el comando `cd` de ChimeraX.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amidas_lib import exige_trayectoria  # noqa: E402

OUT = "results/md/tables/"
s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
exige_trayectoria(s)   # falla antes de sobrescribir nada
lig = s.residues.filter(s.residues.names == "LIG")[0]
lig_heavy = lig.atoms.filter(lig.atoms.element_names != "H")
prot = s.residues.filter(s.residues.polymer_types == 1)
prot_atoms = prot.atoms.filter(prot.atoms.element_names != "H")
ids = s.coordset_ids
ref = None
hb_count, contact_count = Counter(), Counter()
rows = []
for cid in ids:
    s.active_coordset_id = cid
    xyz = lig_heavy.scene_coords
    if ref is None:
        ref = xyz.copy()
    rmsd = float(np.sqrt(((xyz - ref) ** 2).sum(1).mean()))
    hbs = find_hbonds(session, [s], inter_model=False, intra_model=True, cache_da=True)
    partners = set()
    for d, a in hbs:
        if (d.residue is lig) != (a.residue is lig):
            other = a.residue if d.residue is lig else d.residue
            if other in prot:
                partners.add(f"{other.name}{other.number}")
    hb_count.update(partners)
    pxyz = prot_atoms.scene_coords
    dmin = np.full(len(prot_atoms), np.inf)
    for p in xyz:
        dmin = np.minimum(dmin, np.linalg.norm(pxyz - p, axis=1))
    near = {f"{r.name}{r.number}" for r in prot_atoms.filter(dmin < 4.0).unique_residues}
    contact_count.update(near)
    rows.append((len(rows) * 0.05, rmsd, len(partners), ";".join(sorted(partners))))
s.active_coordset_id = ids[0]
n = len(ids)
with open(OUT + "por_fotograma.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tiempo_ns", "rmsd_lig_A", "n_hbonds", "residuos_hbond"]); w.writerows(rows)
with open(OUT + "ocupacion.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["residuo", "hbond_pct", "contacto_pct"])
    # El nombre desempata: sin el, el orden de iteracion del set (que Python
    # aleatoriza) decide entre residuos con el mismo % y el fichero sale
    # distinto en cada ejecucion.
    for r in sorted(set(hb_count) | set(contact_count),
                    key=lambda k: (-contact_count[k], k)):
        w.writerow([r, round(100 * hb_count[r] / n, 1), round(100 * contact_count[r] / n, 1)])
r = np.array([x[1] for x in rows])
session.logger.info(f"frames={n} rmsd_mean={r.mean():.2f} rmsd_max={r.max():.2f} rmsd_last5ns={r[-100:].mean():.2f}")
session.logger.info("HB: " + ", ".join(f"{k} {100*v/n:.0f}%" for k, v in hb_count.most_common(8)))
session.logger.info("CT: " + ", ".join(f"{k} {100*v/n:.0f}%" for k, v in contact_count.most_common(15)))
