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

---

## Flujo de trabajo

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