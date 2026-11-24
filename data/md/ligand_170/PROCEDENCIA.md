# De dónde sale cada fichero de entrada

La simulación del 170 se montó en un directorio aparte y se ha traído aquí
después. Esta tabla dice de dónde salió cada entrada; todas estaban ya en el
repositorio, así que ahora el pipeline las lee directamente de su sitio y no
hay copias congeladas.

| Fichero | Dónde está en el repositorio | Qué es |
|---|---|---|
| `7AAA_repaired.pdb` | `data/md/topology/` | PARP1 reparada con PDBFixer: solo cadena A, sin reconstruir los residuos ausentes (Gly660 y Ser661 quedan fuera a propósito). Es exactamente la misma proteína que usó la simulación del ligando 453, así que las dos son comparables. |
| `docked_ligando_170.mol2` | `results/docking/candidatos/` | Las 10 poses que devolvió GNINA. Se usa la **pose 1**, que es la que seleccionó el cribado. |
| `docked_ligando_170_complex.pdb` | `results/docking/candidatos_complejo/` | El complejo proteína-pose, solo como referencia para comprobar contactos. No entra en la simulación. Ojo: está renumerado, los residuos van con un desplazamiento de −661 (Gly863 aparece como Gly202). |
| el representante estructural | registro **170** de `results/docking/representantes_docking.sdf` | **Sin hidrógenos**. Solo se usa para los órdenes de enlace: las coordenadas que cuentan son las del mol2. El número del ligando es su posición en ese fichero, y `preparar_ligando.py` lo extrae de ahí. |

Lo que sí se guarda en esta carpeta es la salida: `ligando_170.sdf` (la pose
con hidrógenos), `gromacs/` (la topología que se usó en la simulación) y
`frame0.pdb`. `traj_fit.xtc` va en el release, no en el repositorio.

## Qué es el ligando 170

```
Molport-053-317-044
C19H22N6S, 366,5 Da, 26 átomos pesados, 5 enlaces rotables
SMILES: Cc1ncsc1CN1CCC[C@@H](c2ccnc(Nc3cnccn3)c2)C1
```

## Por qué este candidato y no otro

De los 31 candidatos del docking, es **el único que cumple las dos cosas a la
vez**: su pose hace los dos anclajes canónicos del subsitio de nicotinamida y
está poco tensionada.

```
                        tension   rot    Vina   CNN  Intramol  GLY863  SER904
ligando_170                 7,5     5  -10,04  0,54     -0,40    3,41    3,16
ligando_453 (el anterior)  26,9    10  -10,61  0,83     +1,51    3,28    3,31
```

La tensión es la suma, sobre los enlaces rotables, de la energía de torsión de
la pose por encima del mínimo de cada enlace, con OpenFF 2.2.1 — el mismo campo
de fuerza de la simulación. El 453 resultó estar en el 10 % más tensionado de
los 31, y la dinámica le desmontó la pose a los 28 ns.

**Lo que no es bueno de este candidato:** es el peor de los 31 en `Vina_score`
(−10,04, justo en el umbral de −10) y tiene un `CNNscore` bajo (0,54). Por los
criterios originales del cribado es el más flojo que pasa el corte. Eso hace
que la simulación discrimine: si lo que importaba eran las métricas de docking,
esta pose debería fallar; si lo que importa es lo que aprendimos del 453,
debería aguantar.

**El enlace que hay que vigilar** es el CH₂−N de la piperidina, que lleva 5,3
kcal/mol de una barrera de 9,0 (el 59 %). Es el más tensionado de los cinco. Si
la pose se mueve, es el primer sospechoso.

El análisis completo de los 31 candidatos y las limitaciones de estas cifras
están en [docs/limitations.md](../../../docs/limitations.md).
