# Resultados de dinámica molecular

Una carpeta por ligando simulado. Cada una tiene sus figuras, sus tablas
y el log de GROMACS de la producción.

| Ligando | Identificador | Qué pasa en 50 ns | Resultados |
|---|---|---|---|
| **453** | MolPort-000-056-093 | **Cambia de modo de unión a los 27,5 ns.** El RMSD del ligando sube de 1,2 a 5,6 Å y el benzotriazol deja el sándwich de tirosinas | [ligand_453/](ligand_453/README.md) |
| **170** | Molport-053-317-044 | **Se queda en la pose del docking.** RMSD por debajo de 2,4 Å los 50 ns, y a los 36 ns forma los puentes de H de Ser904 y Gly863 | [ligand_170/](ligand_170/README.md) |

`comparacion/` tiene el RMSD del ligando y del backbone del 453
recalculados con el mismo paso de 100 ps que el 170, que es lo que
permite superponer las dos simulaciones en una misma gráfica.

Las dos son **una sola réplica** de 50 ns cada una, trabajo exploratorio
posterior al TFM. La diferencia entre ellas es la de una tirada contra
otra tirada, no una medida: ver [limitations.md](../../docs/limitations.md).

## El protocolo es el mismo

Misma proteína (7AAA cadena A, reparada igual), mismo campo de fuerza
(amber19sb + TIP3P para el sistema, OpenFF 2.2.1 con cargas NAGL para el
ligando), misma caja, mismo equilibrado y misma producción de 50 ns. Los
parámetros están en [02_molecular_dynamics.md](../../docs/02_molecular_dynamics.md).

Las dos simulaciones se lanzaron de manera distinta, y eso sí difiere:

|  | Ligando 453 | Ligando 170 |
|---|---|---|
| Cómo se ejecutó | a mano, paso a paso; la receta se recuperó después del historial de terminal (`scripts/md/run_md.sh`) | con el script unificado `scripts/md/ligando/run_md.sh` |
| Semilla del NVT | `gen-seed = -1` (aleatoria; las efectivas, en [provenance.md](../../docs/provenance.md)) | `gen-seed = 170`, escrita en el `.mdp` |
| Hardware | 4 hilos de CPU, ~10 ns/día | 4 hilos de CPU los primeros 3,47 ns y GPU (CUDA) el resto, 45,6 ns/día |
| Checkpoint | cada 60 min | cada 15 min |
| `.xtc` del equilibrado | no se guardaron | se conservan |
| Fotograma guardado | cada 10 ps | cada 10 ps |
| Trayectoria ligera | 1 cada 50 ps (1001 fotogramas) | 1 cada 100 ps (501 fotogramas) |

Ese último punto es el que hay que recordar al usar los scripts de
`scripts/md_analysis/`: el paso entre fotogramas es su segundo argumento,
y por defecto vale 0,05 ns, que es el del 453.

```bash
runscript scripts/md_analysis/analiza.py          # 453, 50 ps
runscript scripts/md_analysis/analiza.py 170 0.1  # 170, 100 ps
```
