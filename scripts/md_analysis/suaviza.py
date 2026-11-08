# Suaviza la trayectoria cargada en #1 con una media móvil (solo para visualización).
import numpy as np
from chimerax.atomic import all_atomic_structures

VENTANA = 7
s = [m for m in all_atomic_structures(session) if m.id == (1,)][0]
ids = list(s.coordset_ids)
xyz = np.array([s.coordset(i).xyzs for i in ids])
pad = VENTANA // 2
xp = np.concatenate([xyz[:1].repeat(pad, 0), xyz, xyz[-1:].repeat(pad, 0)])
acum = np.cumsum(np.concatenate([np.zeros_like(xp[:1]), xp]), axis=0)
suave = (acum[VENTANA:] - acum[:-VENTANA]) / VENTANA
s.add_coordsets(suave, replace=True)
session.logger.info(f"Trayectoria suavizada: {len(ids)} fotogramas, ventana {VENTANA}")
