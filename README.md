# PARP1 Virtual Screening Pipeline

Pipeline bioinformático para el cribado virtual de inhibidores de PARP1 (Poli ADP-ribosa polimerasa 1) mediante farmacóforos estructurales basados en el bolsillo catalítico y en las interacciones observadas con inhibidores conocidos. Desarrollado como Trabajo de Fin de Máster (TFM) del Máster en Bioinformática.

---

## Descripción del proyecto

PARP1 es una enzima nuclear implicada en la detección y señalización del daño en el ADN. Su dominio catalítico (CAT), compuesto por los subdominios HD (helical domain) y ART (ADP-ribosyl transferase), constituye la principal diana terapéutica de los inhibidores clínicos de PARP (PARPi), como olaparib, niraparib o veliparib. Estos compuestos actúan principalmente mediante la unión competitiva al sitio catalítico de la enzima, interfiriendo con la utilización del cofactor NAD+.

Este proyecto desarrolla un pipeline reproducible y basado principalmente en herramientas de código abierto para realizar un cribado virtual de posibles inhibidores de PARP1.

El objetivo principal no es realizar directamente el descubrimiento de un nuevo fármaco, sino desarrollar y evaluar un flujo de trabajo computacional reproducible capaz de **priorizar candidatos potenciales mediante la integración de información estructural, farmacofórica y de docking molecular**.

El pipeline comprende las siguientes etapas:

1. Recopilación y preparación de estructuras experimentales de PARP1 y de un modelo estructural predicho.
2. Análisis de las interacciones proteína-ligando presentes en estructuras cristalográficas.
3. Generación de un farmacóforo de consenso basado en las interacciones observadas en el bolsillo catalítico.
4. Refinamiento y validación retrospectiva del farmacóforo mediante un conjunto de inhibidores conocidos y moléculas decoy.
5. Cribado virtual de una biblioteca molecular mediante Pharmit.
6. Selección de candidatos y reducción de la redundancia estructural mediante fingerprints moleculares y clustering.
7. Docking molecular de los representantes estructurales seleccionados mediante GNINA.
8. Priorización de candidatos en función de las puntuaciones de docking y de la evaluación mediante la red neuronal de GNINA.
9. Análisis detallado de las interacciones proteína-ligando de los candidatos seleccionados.
10. Comparación de las interacciones observadas en los candidatos con las características definidas por el farmacóforo de consenso.

---

## Flujo de trabajo

El flujo general desarrollado en el proyecto puede resumirse de la siguiente manera:

```text
Estructuras PARP1
       │
       ▼
Preparación de estructuras
       │
       ▼
Análisis de interacciones
       │
       ▼
Farmacóforo de consenso
       │
       ▼
Validación retrospectiva
       │
       ▼
Refinamiento del farmacóforo
       │
       ▼
Cribado virtual con Pharmit
       │
       ▼
1000 moléculas seleccionadas
       │
       ▼
Morgan fingerprints + Butina clustering
       │
       ▼
855 representantes estructurales
       │
       ▼
Docking molecular con GNINA
       │
       ▼
Priorización mediante Vina score + CNNscore
       │
       ▼
Análisis estructural e interacciones
       │
       ▼
Priorización final de candidatos
```

---

## Estructura del repositorio
```
PARP1-Virtual-Screening/
├── info_adicional/
│   └── Información auxiliar del proyecto, estructuras de referencia,
│       datos de validación y documentación del conjunto de estructuras.
│
├── results/
│   ├── FARMACOFORO/
│   │   └── Resultados relacionados con la generación y análisis
│   │       de los farmacóforos y las interacciones proteína-ligando.
│   │
│   ├── PHARMIT/
│   │   └── Resultados del cribado virtual mediante Pharmit,
│   │       incluyendo los candidatos seleccionados y los datos
│   │       utilizados para su posterior análisis.
│   │
│   └── DOCKING/
│       └── Resultados del docking molecular realizado mediante GNINA,
│           incluyendo las estructuras y datos empleados para
│           la selección y análisis de candidatos.
│
├── scripts/
│   └── Scripts Python que automatizan las diferentes etapas del pipeline,
│       desde la preparación de estructuras hasta el clustering y docking.
│
├── structures/
│   ├── raw/
│   │   └── Estructuras PDB originales obtenidas de las bases de datos.
│   │
│   ├── processed/
│   │   └── Estructuras procesadas y preparadas para los análisis.
│   │
│   └── contacts/
│       └── Resultados del análisis automatizado de las interacciones
│           proteína-ligando mediante PLIP.
│
├── parp1_pipeline.yml
│   └── Definición del entorno Conda utilizado en el proyecto.
│
└── README.md
    └── Documentación del proyecto.
```
> **Nota:** El directorio `redocking/` se utiliza para la validación del protocolo de docking mediante el redocking de estructuras cristalográficas. Debido al elevado volumen de datos generado durante esta etapa, este directorio no se incluye en la versión final del repositorio de GitHub.

> **Nota:** Los archivos estructurales de gran tamaño pueden estar excluidos del control de versiones mediante .gitignore. Las diferentes etapas del pipeline pueden reproducirse mediante los scripts y los datos auxiliares proporcionados.

---

## Estructras utilizadas

Todas las estructuras experimentales utilizadas corresponden al **dominio catalítico (CAT) de PARP1 humana**, compuesto por los subdominios HD (`661–786`) y ART (`786–1014`), con la excepción del modelo de longitud completa generado mediante AlphaFold.

La información relativa a las estructuras PDB empleadas, sus ligandos y los datos experimentales asociados puede consultarse en la [tabla de estructuras PDB](./info_adicional/tabla_pdb.csv).

El modelo de longitud completa de PARP1 se obtuvo de la [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/entry/AF-P09874-F1) y se utilizó como referencia estructural adicional.

Las estructuras correspondientes a los PDB 9ETQ y 9ETR no se encuentran disponibles para su descarga en formato PDB convencional.

---

## Análisis de interacciones proteína-ligando

Las interacciones entre los inhibidores conocidos y el bolsillo catalítico de PARP1 se analizaron mediante [PLIP (*Protein-Ligand Interaction Profiler*)](https://doi.org/10.1093/nar/gkaf361)

El análisis permitió identifciar los residuos y tipos de interacción más frecuentes entre los inhibidores y la proteína, proporcionando la base estructural utilizada posteriormente para la construcción del farmacóforo consenso.

Entre las interaciones identificadas se encuentran:

* Puentes de hidrógeno (tanto dondadores como aceptores)
* Interacciones hidrofóbicas
* Interacciones π–π
* Otras interacciones que se han omitido en este trabajo debido a que no son relevantes (interacciones iónicas)

Los resultados estructurados de las interaccciones pueden consultarse en [results/FARMACOFORO/interacciones.csv](./results/FARMACOFORO/interacciones.csv)

---

## Construcción del farmacóforo

A partir de la interacciones identificadas en las estructuras experimentales inlcuidas en el [*Training set*](./info_adicional/tabla_pdb.csv), se construyeron modelos faramcofóricos que representan las características estructurales relevante para la interación con el bolsillo catalítico de PARP1.

Las características farmacofóricas se definieron a partir de la frecuencia y conservación de las interacciones observadas entre los diferente complejos proteína-ligando.

Se generaron así diferentes niveles de exigencia:

* **Core:** características consideradas esenciales
* **Recommended:** características altamente conservadas pero no estrictamente obligatorias
* **All:** conjunto completo de las características identificadas en todas las estructuras

Los modelos farmacoforicos y los datos empleados se encuentra en: [results/FARMACOFORO](./results/FARMACOFORO/)

---

## Validación retrospectiva y refinamiento del farmacóforo

Los modelos farmacofóricos iniciales fueron evaluados mediante una **validación retrospectiva**, utilizando un conjunto formado por inhibidores conocidos y moléculas decoy.

Esta etapa permitió evaluar la capacidad del farmacóforo para discriminar entre moléculas activas y moléculas que no presentan las características estructurales esperadas.

A partir de los resultados obtenidos se realizó un refinamiento de las características farmacofóricas, seleccionando aquellas que proporcionaban una representación más adecuada de los patrones de reconocimiento observados experimentalmente.

Los datos utilizados para esta etapa y los modelos refinados se encuentran en: [results/PHARMIT/farmacoforo_refined](./results/PHARMIT/farmacoforo_refined.csv)

---

## Cribado virtual mediante Pharmit

Los farmacóforos refinados se utilizaron para realizar un cribado virtual mediante [Pharmit](https://doi.org/10.1093/nar/gkw287).

El cribado permitió identificar moléculas compatibles con las características farmacofóricas definidas para el bolsillo catalítico de PARP1.

A partir de los resultados obtenidos se seleccionaron inicialmente las **1000 moléculas con mejores resultados de ajuste al farmacóforo (RMSD)**, que fueron utilizadas como conjunto de partida para la reducción de la redundancia estructural.

Los resultados del cribado se encuentran en: [results/PHARMIT/](./results/FARMACOFORO)

---

## Clustering y selección de representantes estructurales

Las 1000 moléculas seleccionadas mediante Pharmit fueron sometidas a un proceso de clustering molecular con el objetivo de reducir la redundancia estructural antes del docking.

Para representar las moléculas se utilizaron [**fingerprints de Morgan**](https://doi.org/10.1021/ci100050t), con:

* Radio: `2`
* Tamaño: `2048 bits`

La similitud estructural entre moléculas se evaluó mediante el coeficiente de Tanimoto y posteriormente se realizó clustering utilizando el [**algoritmo de Butina**](https://doi.org/10.1021/ci9803381).

Se utilizó una distancia de clustering de `0.2`, equivalente a considerar una similitud de Tanimoto de aproximadamente `0.8`.

Este procedimiento permitió reducir el conjunto inicial de 1000 moléculas a **855 representantes estructurales**, que fueron utilizados posteriormente como conjunto de entrada para el docking molecular.

El proceso de selección y clustering se automatiza mediante: [scripts/10_clustering_repr.py](./scripts/10_clustering_repr.py)

---

## Docking molecular

Las **855 moléculas representativas** obtenidas tras el proceso de clustering se sometieron a docking molecular mediante [GNINA](https://doi.org/10.1186/s13321-025-00973-x).

GNINA combina métodos de docking clásicos derivados de AutoDock Vina con modelos de aprendizaje profundo para evaluar las poses generadas.

Para cada ligando se generaron múltiples poses de unión y se obtuvieron diferentes métricas:

* *`Vina_score`*: estimación energética de la afinidad de unión de la pose.
* *`CNNscore`*: estimación de la plausibilidad de la pose mediante el modelo de aprendizaje profundo.
* *`CNNaffinity`*: estimación adicional de la afinidad mediante el modelo CNN.

El docking se automatiza mediante: [scripts/11_docking.py](./scripts/11_docking.py)

Los resultados generados se encuentran en: [results/DOCKING/](./results/DOCKING/)

---