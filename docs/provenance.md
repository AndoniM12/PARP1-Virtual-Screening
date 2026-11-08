# Procedencia y lagunas

Qué se ha podido recuperar de los ficheros del proyecto, qué se ha
recuperado del historial de terminal y qué no se ha recuperado en
absoluto.

Nada de lo que aparece como «no recuperado» se ha reconstruido de
memoria ni inventado.

## Recuperado del historial de terminal

Estos pasos **no existían como script**. Solo estaban en el historial
de la terminal y se han convertido en código:

| Paso | Dónde está ahora |
|---|---|
| Montaje del complejo proteína-ligando (concatenar los dos `.gro`) | `scripts/md/montar_complejo.py` |
| Inserción de `LIG 1` en `[ molecules ]` de la topología | `scripts/md/run_md.sh`, paso 8 |
| `gmx editconf -bt dodecahedron -d 1.0` | `run_md.sh`, paso 9 |
| `gmx solvate -cs spc216.gro` | ídem |
| `gmx grompp ... -maxwarn 1` y `gmx genion -neutral` | ídem |
| Los `grompp` de minimización, NVT, NPT y producción | `run_md.sh`, pasos 10-13 |
| `gmx trjconv -pbc mol -center` | `run_md.sh`, paso 14 |

El script original de montaje del complejo llevaba `assert` con los
números de átomos fijados (5546 y 53). Se han conservado como
comprobación, pero avisando en lugar de abortar.

## No recuperado

### 1. El `trjconv` que generó `traj_ligera.xtc`

**Lo más importante de esta lista**, porque de ahí depende cómo se lee
la gráfica de RMSD.

Lo que se sabe, del README del análisis original:

- sin agua
- ajustada al backbone
- 1 fotograma cada 50 ps, 1001 en total (50 ns)

Lo que **no** se sabe: el comando exacto. En concreto, sobre qué grupo
del índice se hizo el ajuste y con qué opción de `-fit`
(`rot+trans`, `progressive`…). El `index.ndx` de la ejecución original
sí existe (`simulation/50ns/index.ndx`, 1,5 MB), pero no consta qué
grupo se eligió.

El historial sí registra `gmx make_ndx -f md_50ns.tpr -o index.ndx` y
`gmx trjconv ... -pbc mol -center`, pero ese último produjo
`md_50ns_centered.xtc`, que no es `traj_ligera.xtc`: el centrado no
incluye `-fit` ni quita el agua ni reduce a 50 ps.

**Consecuencia:** la gráfica de RMSD se puede volver a generar a partir
de `traj_ligera.xtc` (publicada como release), pero no se puede
reconstruir `traj_ligera.xtc` desde cero de forma idéntica sin volver a
decidir ese comando.

### 2. El binario de GNINA

La **versión sí consta**: GNINA **1.3.3**, documentada en la tabla de
herramientas de [01_virtual_screening.md](01_virtual_screening.md).

Lo que ya no está es el binario: GNINA se instaló al margen de conda y
no queda en el sistema (el entorno `gnina` existe pero está vacío de
él). Tampoco hay ningún log de GNINA versionado que registre la versión
del modelo de red neuronal usado para el CNN score.

Los parámetros sí están en `scripts/screening/11_docking.py`:

- caja: centro (8,59, 39,35, 9,53), tamaño 22 × 16 × 24 Å
- `--exhaustiveness 8`, `--num_modes 10`, **`--seed 1`**

La semilla está fijada y la versión del programa se conoce, así que el
docking es reproducible salvo por la versión del modelo de red neuronal
usado para el CNN score, que GNINA elige por defecto y no quedó
registrada.

### 3. La biblioteca cribada con Pharmit

El cribado se hizo en la interfaz web de Pharmit. No quedó registrada
ni la versión de la base de datos ni la fecha de la consulta.

La salida completa (`pharmit.sdf`, **2,8 GB, 1 479 579 moléculas**) **no
está en el repositorio** por su tamaño. Estuvo a punto de perderse: no
quedaba ninguna copia en disco y solo sobrevivía como un objeto suelto,
nunca confirmado, dentro de `.git` del repositorio original. Se ha
extraído a `~/Documents/pharmit_backup/pharmit.sdf`.

Lo que sí está versionado es el resultado filtrado:
`results/pharmit/top1000_rmsd.sdf`.

### 4. Las semillas de la dinámica

Los `.mdp` usan `gen-seed = -1`, es decir, semilla aleatoria. La
simulación original **no es reproducible bit a bit** salvo que se fijen.

Las semillas efectivas sí constan en las cabeceras de los `.log`:

| Etapa | `ld-seed` |
|---|---|
| Minimización | 2078777343 |
| NVT | −553681930 |
| NPT | 1475263481 |
| Producción 50 ns | 732815095 |

Para repetir la ejecución original hay que ponerlas a mano en los
`.mdp` correspondientes.

### 5. La versión de ChimeraX

Deducida del nombre del paquete instalado que figura en el historial
(`ucsf-chimerax_1.12ubuntu24.04_amd64.deb`) → **1.12**. No estaba
declarada en ningún fichero del proyecto.

## PLIP no es determinista

Comprobado al reejecutar el pipeline desde cero partiendo solo de lo publicado.

PLIP da resultados distintos entre ejecuciones sobre un fichero de entrada
**idéntico**. Diez ejecuciones sobre `9ETQ_chainA.pdb`:

```
3 4 3 4 3 4 3 3 4 3     enlaces de hidrógeno detectados
- S - S - S - - S -     ¿aparece Arg878?
```

El puente de hidrógeno con **Arg878 aparece en 4 de cada 10 ejecuciones**. Es un
contacto en el límite geométrico. La ejecución publicada lo capturó; la mayoría
no lo hace.

Lo que arrastra, cuando no se detecta:

| ARG878 HBA | PDBs | x, y, z | Radio | Conservación |
|---|---|---|---|---|
| Publicado | 2 | 2,519 / 42,647 / 10,828 | 3,10 | 0,33 |
| Al reejecutar | 1 | 4,624 / 43,117 / 10,436 | 1,00 | 0,17 |

En el modelo `all` refinado ese punto se desplaza unos **2,2 Å** y su radio baja
de 1,2 a 1,0.

### Alcance real: no afecta al cribado

Conviene acotarlo, porque suena peor de lo que es.

Arg878 HBA está clasificada como **Optional**, así que aparece **únicamente en
el modelo `all`**. Y `all` tiene una tasa de recuperación del **0 %** en la
validación retrospectiva (`results/pharmit/resultados_pharmit.csv`): no recupera
ninguno de los inhibidores conocidos, de modo que no sirve para cribar y no se
usó para ello.

Los modelos que sí se usaron son los que recuperan activos:

| Modelo | Características | Recuperación (refinado) | ¿Contiene Arg878? | ¿Reproduce? |
|---|---|---|---|---|
| `core` | 3 (Gly863 HBD, Gly863 HBA, Tyr907 aromática) | **100 %** | No | **Idéntico** |
| `recommended` | 4 (las anteriores + Tyr889 hidrofóbica) | 82,61 % | No | **Idéntico** |
| `all` | 12 | 0 % | Sí | Difiere en ese punto |

Es decir: **la no-determinación de PLIP no toca ninguna característica de los
farmacóforos con los que se cribó**. El cribado virtual, los 855 representantes
y el docking son reproducibles sin reservas.

Queda documentado porque es un comportamiento real de la herramienta que
conviene conocer, no porque comprometa los resultados.

**No se han regenerado los resultados publicados con una ejecución nueva**: los
que están en el repositorio son los de la ejecución original, que es el registro
real del trabajo. Quien reproduzca el pipeline obtendrá casi siempre la variante
sin Arg878.

### Qué sí reproduce exacto

Verificado partiendo de un `git archive` del repositorio:

| Paso | Resultado |
|---|---|
| 02, estructuras procesadas | **las 31 byte a byte idénticas** |
| 03, contactos | 6 estructuras |
| 05, tabla de consenso | 16 características |
| 06, modelos sin refinar | `core` y `recommended` idénticos; `all` difiere solo en Arg878 |
| 09, modelos refinados | `core` (3) y `recommended` (4) idénticos; `all` (12) difiere solo en Arg878 |
| 10, representantes | **855, mismos nombres y mismos SMILES** |
| `reparar_proteina.py` (fase 2) | idéntica átomo por átomo; solo cambia la fecha de un `REMARK` |

El docking no se reejecutó: con los 855 ligandos de entrada idénticos, GNINA
daría el mismo resultado.

## Un error corregido al reorganizar

`07_export_pharmit.py` construía los tres modelos de farmacóforo filtrando
`farmacoforo_refined.csv` por la columna `Classification`. Eso solo es válido si
la tabla es por modelo, y no lo es: es la **unión** de los tres.

El resultado era un modelo `all` con **13** características en vez de **12**,
que colaba la hidrofóbica de Tyr889 — la que el refinamiento elimina
precisamente porque la aromática de Tyr896 ya la cubre.

Los JSON que estaban guardados en el repositorio sí eran correctos (12/3/4),
pero eran **anteriores** a los dos scripts de refinamiento (11:19–11:54 frente a
12:36 y 13:30 del 8 de octubre de 2026), así que el pipeline ya no sabía
regenerarlos.

La exportación vive ahora dentro de `09_refinamiento_farmacoforo.py`, que la hace desde la
composición por modelo, que es la única correcta. Se ha comprobado que los JSON
que genera son idénticos, coordenada a coordenada, a los que estaban guardados.
`07_export_pharmit.py` se ha eliminado; sigue recuperable del historial de git.

La variante anterior se conserva como `09b_refinamiento_por_modelo.py`: el
mismo refinamiento, con una salida distinta (una fila por característica y
modelo, con columna `Model`) que sirve para inspeccionar. Escribe a
`farmacoforo_refined_por_modelo.csv` para no pisar la salida del script vigente.

Los dos scripts se llamaban antes `09_refinamiento_farmacoforo.py` y
`refinamiento.py`; este último era el vigente pese al nombre, así que se han
intercambiado para que el numerado sea el del pipeline.

## Rutas que se han cambiado

Los scripts originales llevaban rutas absolutas y rutas relativas a una
estructura de carpetas que ya no existe. Se han adaptado para poder
lanzarlos desde la raíz del repositorio:

| Antes | Ahora |
|---|---|
| `/home/andonim12/Documents/PARP1/analisis_ligando453_chimerax/resultados/tablas/` | `results/md/tables/` |
| `/home/andonim12/Documents/PARP1/analisis_ligando453_chimerax/datos/` | `data/md/` |
| `/home/andonim12/Documents/PARP1/analisis_ligando453_chimerax/render/` | `run/render/` |
| `/home/andonim12/miniconda3/envs/openff/.../openff-gnn-am1bcc-1.0.0.pt` | se resuelve desde el paquete instalado |
| `$HOME/Documents/PARP1-Molecular-Dynamics/simulation/50ns/` | `run/prod/` |
| `preparation/ligand/ligand_453/` | `data/md/ligand_453/` (entrada) y `run/ligand/` (salida) |
| `data/7AAA.pdb` | `data/raw/7AAA.pdb` |
| `data/p-l_complex_pdb/` | `results/docking/candidatos_complejo/` |
| `structures/{raw,processed,contacts}` | `data/structures/{raw,processed,contacts}` |
| `info_adicional/` | `data/raw/` |
| `results/FARMACOFORO/`, `results/farmacoforo.csv` | `results/pharmacophore/` |
| `results/PHARMIT/` | `results/pharmit/` |
| `results/DOCKING/` | `results/docking/` |
| `results/interacciones.csv` | `results/pharmacophore/interacciones.csv` |

`data/md/topology/topol.top` conserva sus directivas `#include` tal como
estaban en la ejecución original: apuntan a la estructura de carpetas
antigua. **No se ha modificado a propósito**, por ser evidencia de lo
que se ejecutó. Al reproducir la simulación, `run_md.sh` regenera la
topología desde cero con `pdb2gmx`.

## Ficheros originales

Nada de esto se ha movido ni borrado. Los directorios originales siguen
intactos en la máquina donde se hizo el trabajo:

- `PARP1-Virtual-Screening/` — repositorio del TFM
- `Molecular-Dynamics_(prueba)/` — dinámica con GROMACS, incluidas las
  trayectorias completas (1,8 GB) que no se han traído aquí
- `analisis_ligando453_chimerax/` — análisis, vídeos y
  `traj_ligera.xtc`

Este repositorio es una copia reorganizada, no un traslado.
