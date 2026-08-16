# PARP1 Virtual Screening Pipeline

Pipeline bioinformático para el cribado virtual de inhibidores de PARP1 (Poli ADP-ribosa polimerasa 1) mediante farmacóforos estructurales basados en el bolsillo catalítico y en las interacciones observadas con inhibidores conocidos. Desarrollado como Trabajo de Fin de Máster (TFM).

---

## Descripción del proyecto

PARP1 es una enzima nuclear implicada en la detección y señalización del daño en el ADN. Su dominio catalítico (CAT), compuesto por los subdominios HD y ART, es la principal diana terapéutica de los inhibidores clínicos de PARP (PARPi), como olaparib, niraparib o veliparib, cuyo mecanismo de acción se basa en la inhibición competitiva frente al cofactor NAD+.

Este proyecto desarrolla un pipeline reproducible y basado principalmente en herramientas de código abierto para realizar un cribado virtual de posibles inhibidores de PARP1.

El pipeline comprende las siguientes etapas:

1. Recopilación y preparación de estructuras experimentales de PARP1 y de un modelo estructural predicho.
2. Análisis de las interacciones proteína-ligando presentes en estructuras cristalográficas.
3. Generación de un farmacóforo de consenso basado en las interacciones observadas en el bolsillo catalítico.
4. Refinamiento y validación retrospectiva del farmacóforo mediante un conjunto de inhibidores conocidos y moléculas decoy.
5. Cribado virtual de una biblioteca molecular mediante Pharmit.
6. Selección de candidatos y reducción de la redundancia estructural mediante fingerprints moleculares y clustering.
7. Docking molecular de los representantes estructurales seleccionados mediante GNINA.
8. Priorización de los candidatos en función de las puntuaciones de docking y de la evaluación mediante la red neuronal de GNINA.
9. Análisis detallado de los candidatos mejor puntuados mediante el estudio de sus interacciones proteína-ligando.

El objetivo del proyecto es desarrollar y evaluar un flujo de trabajo reproducible para la priorización computacional de candidatos, y no realizar por sí mismo el desarrollo de un nuevo fármaco.

---

## Estructura del repositorio

```text
PARP1-Virtual-Screening/
├── info_adicional/
│   └── Información auxiliar del proyecto, estructuras de referencia,
│       datos de validación y documentación del conjunto de estructuras.
│
├── results/
│   ├── FARMACOFORO/
│   │   └── Resultados relacionados con la generación y análisis
│   │       de los farmacóforos.
│   │
│   ├── PHARMIT/
│   │   └── Resultados del cribado virtual mediante Pharmit,
│   │       incluyendo los candidatos seleccionados y los datos
│   │       utilizados para su posterior análisis.
│   │
│   └── DOCKING/
│       └── Resultados del docking molecular realizado mediante GNINA,
│           incluyendo estructuras, resultados y datos empleados
│           para la selección de candidatos.
│
├── scripts/
│   └── Scripts Python que automatizan las diferentes etapas del pipeline,
│       desde la preparación de estructuras hasta el clustering y docking.
│
├── structures/
│   ├── raw/
│   │   └── Estructuras PDB originales.
│   │
│   ├── processed/
│   │   └── Estructuras procesadas y preparadas para los análisis.
│   │
│   └── contacts/
│       └── Resultados del análisis automatizado de las interacciones
│           proteína-ligando.
│
├── parp1_pipeline.yml
│   └── Definición del entorno Conda utilizado.
│
└── README.md
    └── Documentación del proyecto.
```

> **Nota:** Las estructuras y otros archivos de gran tamaño pueden estar excluidos del control de versiones mediante `.gitignore`. Las diferentes etapas del pipeline pueden reproducirse mediante los scripts y los datos auxiliares proporcionados.

> **Nota:** El directorio `redocking/` se utiliza para la validación mediante redocking de estructuras experimentales y no se incluye en el repositorio final de GitHub debido al volumen de datos generado.

---

## Estructuras utilizadas

Todas las estructuras experimentales utilizadas corresponden al **dominio catalítico (CAT) de PARP1 humana**, compuesto por los subdominios HD (661–786) y ART (786–1014), con la excepción del modelo de longitud completa generado mediante AlphaFold.

La información relativa a las estructuras PDB empleadas, sus ligandos y los datos experimentales asociados puede consultarse en la [tabla de estructuras PDB](info_adicional/tabla_pdb.csv).

El modelo de longitud completa de PARP1 se obtuvo de [AlphaFold Database](https://alphafold.ebi.ac.uk/) y se utilizó como referencia estructural adicional.

Las estructuras correspondientes a los PDB 9ETQ y 9ETR no se encuentran disponibles para su descarga en formato PDB convencional.

---

## Análisis de interacciones proteína-ligando

Las interacciones entre los inhibidores conocidos y el bolsillo catalítico de PARP1 se analizaron mediante [PLIP (Protein-Ligand Interaction Profiler)](https://plip-tool.biotec.tu-dresden.de/).

El análisis permitió identificar los residuos y tipos de interacción más frecuentes entre los inhibidores y la proteína, proporcionando la base estructural utilizada posteriormente para la construcción del farmacóforo de consenso.

Los resultados estructurados de las interacciones pueden consultarse en [`results/interacciones.csv`](results/interacciones.csv).

---

## Construcción y validación del farmacóforo

A partir de las interacciones identificadas en las estructuras experimentales se construyeron modelos farmacofóricos que representan las características estructurales relevantes para la interacción con el bolsillo catalítico de PARP1.

Se generaron diferentes niveles de exigencia del farmacóforo:

* `Core`
* `Recommended`
* `All`

Los modelos fueron posteriormente evaluados mediante una validación retrospectiva utilizando un conjunto de inhibidores conocidos y moléculas decoy.

Los resultados y modelos generados durante esta etapa se encuentran en el directorio [`results/FARMACOFORO/`](results/FARMACOFORO/).

---

## Cribado virtual mediante Pharmit

Los farmacóforos refinados se utilizaron para realizar un cribado virtual mediante [Pharmit](https://pharmit.csb.pitt.edu/).

El cribado permitió seleccionar inicialmente un conjunto de moléculas compatibles con las características farmacofóricas definidas para el bolsillo catalítico de PARP1.

Los resultados del cribado y los datos empleados en las etapas posteriores se encuentran en [`results/PHARMIT/`](results/PHARMIT/).

A partir de las **1000 moléculas seleccionadas inicialmente mediante Pharmit**, se realizó una reducción de la redundancia estructural mediante fingerprints de Morgan y clustering utilizando el algoritmo de Butina.

Este procedimiento permitió obtener **855 moléculas representativas**, que fueron utilizadas posteriormente como conjunto de entrada para el docking molecular.

El proceso de selección y clustering se automatiza mediante el script [`scripts/10_clustering_repr.py`](scripts/10_clustering_repr.py).

---

## Docking molecular

Las **855 moléculas obtenidas tras el proceso de clustering** se sometieron a docking molecular mediante [GNINA](https://github.com/gnina/gnina).

GNINA permite combinar métodos de docking clásicos con modelos de aprendizaje profundo para evaluar las poses generadas. Para cada ligando se generaron múltiples poses y se obtuvieron diferentes métricas:

* `Vina_score`: estimación de la energía de unión de la pose.
* `CNNscore`: estimación de la plausibilidad de la pose mediante el modelo de aprendizaje profundo.
* `CNNaffinity`: estimación adicional de la afinidad mediante el modelo CNN.

El docking se automatiza mediante [`scripts/11_docking.py`](scripts/11_docking.py).

Los resultados generados se encuentran en [`results/DOCKING/`](results/DOCKING/).

La distribución de las puntuaciones obtenidas y la relación entre las diferentes métricas se utilizaron para realizar una primera priorización de los candidatos.

---

## Priorización de candidatos

La selección de candidatos no se basó exclusivamente en una única métrica.

En particular, se analizaron conjuntamente `Vina_score` y `CNNscore`, ya que ambas métricas proporcionan información diferente sobre las poses generadas por GNINA.

El `Vina_score` permite priorizar poses con una estimación energética favorable, mientras que `CNNscore` aporta información adicional sobre la plausibilidad de la pose según el modelo de aprendizaje profundo.

A partir de este análisis se seleccionaron los candidatos con valores simultáneamente favorables de ambas métricas para su posterior estudio estructural.

Los resultados completos del docking y los datos utilizados para la priorización se encuentran en [`results/DOCKING/`](results/DOCKING/).

---

## Redocking y validación del docking

La capacidad del protocolo de docking para reproducir las conformaciones experimentales se evaluó mediante un procedimiento de **redocking**.

Para ello, se utilizaron ligandos procedentes de estructuras cristalográficas de PARP1 y se compararon las poses predichas por GNINA con las conformaciones experimentales mediante RMSD.

Esta etapa permite evaluar la capacidad del protocolo empleado para reproducir las orientaciones observadas experimentalmente antes de aplicarlo al conjunto de candidatos.

Los resultados del redocking se utilizan como parte de la validación del pipeline.

El directorio `redocking/` no se incluye en el repositorio final de GitHub debido al tamaño de los archivos generados.

---

## Análisis estructural de los candidatos

Tras la primera priorización basada en las métricas de GNINA, los candidatos mejor puntuados serán sometidos a un análisis estructural más detallado.

Este análisis permitirá estudiar:

* Interacciones con residuos clave del bolsillo catalítico.
* Formación de enlaces de hidrógeno.
* Interacciones hidrofóbicas.
* Interacciones π–π y catión–π cuando estén presentes.
* Conservación de las características farmacofóricas.
* Similitud de las interacciones con las observadas en inhibidores conocidos.

El objetivo de esta etapa es complementar las puntuaciones obtenidas mediante docking con una interpretación estructural de las poses, permitiendo una priorización final más informada de los candidatos.

---

## Instalación y uso

### 1. Clonar el repositorio

```bash
git clone https://github.com/AndoniM12/PARP1-Virtual-Screening.git
cd PARP1-Virtual-Screening
```

### 2. Crear el entorno Conda

```bash
conda env create -f parp1_pipeline.yml -n parp1_pipeline
conda activate parp1_pipeline
```

### 3. Descargar las estructuras

Las estructuras cristalográficas de PARP1 se obtienen desde [RCSB PDB](https://www.rcsb.org/) y el modelo de longitud completa se obtiene de AlphaFold Database.

```bash
python scripts/01_descargar_pdbs.py
```

Las estructuras originales se almacenan en `structures/raw/`.

### 4. Procesar las estructuras

```bash
python scripts/02_procesar_pdb.py
```

Este paso permite obtener las estructuras procesadas y extraer la cadena utilizada posteriormente en los análisis.

### 5. Identificar las interacciones proteína-ligando

```bash
python scripts/03_contactos.py
```

El análisis se realiza mediante PLIP y genera los resultados de las interacciones en `structures/contacts/`.

### 6. Procesar las interacciones

```bash
python scripts/04_lectura.py
```

Este script recopila los resultados obtenidos mediante PLIP y genera los datos estructurados de las interacciones detectadas.

### 7. Generar los datos del farmacóforo

```bash
python scripts/05_farmacoforo_data.py
```

Este paso genera los datos utilizados para construir el consenso de interacciones.

### 8. Generar los farmacóforos

```bash
python scripts/06_gen_farmacoforo.py
```

Se generan los diferentes modelos farmacofóricos utilizados posteriormente en la validación y el cribado virtual.

### 9. Exportar los farmacóforos a Pharmit

```bash
python scripts/07_export_pharmit.py
```

Permite preparar los modelos farmacofóricos para su utilización en Pharmit.

### 10. Generar el conjunto de validación

```bash
python scripts/08_validation_set.py
```

Prepara el conjunto de moléculas utilizado para la validación retrospectiva del farmacóforo.

### 11. Refinar el farmacóforo

```bash
python scripts/09_refinamiento_farmacoforo.py
```

Permite evaluar y refinar los modelos farmacofóricos utilizando el conjunto de validación.

### 12. Seleccionar representantes estructurales

```bash
python scripts/10_clustering_repr.py
```

Las moléculas seleccionadas mediante Pharmit se representan mediante fingerprints de Morgan y se agrupan mediante clustering de Butina para reducir la redundancia estructural.

### 13. Realizar el docking molecular

```bash
python scripts/11_docking.py
```

Los representantes estructurales seleccionados se someten a docking molecular mediante GNINA.

---

## Entorno de trabajo

| Herramienta  | Versión   | Uso                                                           |
| ------------ | --------- | ------------------------------------------------------------- |
| Python       | 3.11      | Scripts del pipeline                                          |
| Biopython    | 1.87      | Lectura y procesado de estructuras PDB                        |
| NumPy        | 2.4.6     | Cálculos numéricos                                            |
| RDKit        | 2025.03.6 | Procesamiento molecular, fingerprints y clustering            |
| PyMOL        | 3.1.0     | Modificación de estructuras y visualización                   |
| PLIP         | 3.0.0     | Identificación automatizada de interacciones proteína-ligando |
| UCSF Chimera | —         | Inspección y validación visual                                |
| Pharmit      | —         | Construcción del farmacóforo y cribado virtual                |
| GNINA        | 1.3.3     | Docking molecular y evaluación mediante CNN                   |
| AlphaFold DB | v6        | Modelo de longitud completa de PARP1                          |
| RCSB PDB     | —         | Estructuras cristalográficas                                  |

---

## Estado actual del proyecto

* [x] Recopilación de estructuras de PARP1
* [x] Descarga automatizada de estructuras
* [x] Extracción y procesamiento de las estructuras
* [x] Análisis de contactos proteína-ligando mediante PLIP
* [x] Generación del consenso de interacciones
* [x] Construcción de los farmacóforos
* [x] Validación retrospectiva del farmacóforo
* [x] Refinamiento de los farmacóforos
* [x] Cribado virtual mediante Pharmit
* [x] Selección de las 1000 moléculas con mejores resultados
* [x] Clustering mediante fingerprints de Morgan y Butina
* [x] Selección de 855 representantes estructurales
* [x] Docking molecular mediante GNINA
* [x] Análisis inicial de las métricas de docking
* [ ] Análisis detallado de los candidatos mejor puntuados
* [ ] Análisis de interacciones proteína-ligando de los candidatos finales
* [ ] Priorización final de candidatos
* [ ] Documentación final del pipeline

---

## Referencias

* [RCSB Protein Data Bank](https://www.rcsb.org/)
* [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
* [PLIP — Protein-Ligand Interaction Profiler](https://plip-tool.biotec.tu-dresden.de/)
* [Pharmit](https://pharmit.csb.pitt.edu/)
* [GNINA](https://github.com/gnina/gnina)
* [RDKit](https://www.rdkit.org/)
* [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/)
* [PyMOL](https://pymol.org/)

---

## Autor

**Andoni Moreno Lanceta**
Trabajo de Fin de Máster — Bioinformática
