# Renderiza la animación a 60 fps con cámara lenta ×6 en la transición (27-31 ns).
# Los fotogramas intermedios se interpolan linealmente entre coordsets consecutivos.
import csv
import math
import os

import numpy as np
from chimerax.atomic import all_atomic_structures
from chimerax.core.commands import run

BASE = "run/render/"            # en .gitignore
os.makedirs(BASE + "mol", exist_ok=True)
RAPIDO, LENTO = 1.5, 0.25            # coordsets avanzados por fotograma de vídeo
INICIO_LENTO, FIN_LENTO = 541, 621   # 27.0 - 31.0 ns
TEMPORAL = 9999

s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
ids = [c for c in s.coordset_ids if c != TEMPORAL]
ultimo = max(ids)
xyz = {c: s.coordset(c).xyzs.copy() for c in ids}

plan, f = [], 1.0
while f <= ultimo:
    lento = INICIO_LENTO <= f < FIN_LENTO
    plan.append((f, lento))
    f += LENTO if lento else RAPIDO

run(session, "2dlabels delete")  # los textos se añaden después con componer.py

angulo_prev = 0.0
with open("results/md/tables/plan_render_60fps.csv", "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["salida", "coordset", "tiempo_ns", "lento"])
    for i, (f, lento) in enumerate(plan):
        a = int(math.floor(f))
        b = min(a + 1, ultimo)
        x = f - a
        s.add_coordset(TEMPORAL, (1 - x) * xyz[a] + x * xyz[b])
        s.active_coordset_id = TEMPORAL
        angulo = 12 * math.sin(2 * math.pi * i / 600)
        run(session, f"turn y {angulo - angulo_prev:.4f} center #1:LIG", log=False)
        angulo_prev = angulo
        t = (f - 1) * 0.05
        run(session, f"save {BASE}mol/{i:04d}.png width 1600 height 1000 supersample 2", log=False)
        w.writerow([i, f"{f:.2f}", f"{t:.3f}", int(lento)])
        out.flush()
s.active_coordset_id = 1
session.logger.info(f"Renderizados {len(plan)} fotogramas")
