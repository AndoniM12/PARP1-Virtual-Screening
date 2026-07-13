# PARP1 Virtual Screening Pipeline

Pipeline bioinformático para el cribado virtual de inhibidores de PARP1 (Poli ADP-ribosa polimerasa 1) mediante farmacóforos estructurale basado en el bolsillo catálitico y fármacos ya clínicamente probados. Desarrollado como Trabajo de Fin de Máster (TFM).

---

## Descripción del proyecto

PARP1 es una enzima nuclear implicada en la detección y señalización del daño en el ADN. Su dominio catalítico (CAT), compuesto por los subdominios HD y ART, es la diana terapéutica principal de los inhibidores clínicos de PARP (PARPi), como olaparib, niraparib o veliparib, cuyo mecanismo de acción se basa en la inhibición competitiva frente al cofactor NAD+.

Este proyecto construye un pipeline reproducible y de código abierto para:

1. Descargar y preparar estructuras cristalográficas y modelos predichos de PARP1
2. Analizar las interacciones proteína-ligando en el bolsillo catalítico
3. Construir un farmacóforo de consenso basado en estructuras ya publicadas y en el bolsillo catalítico
4. Realizar un cribado virtual contra librerías de moléculas comerciales (ZINC o PubCHEM)
5. Validar el farmacóforo mediante validación retrospectiva con activos conocidos y decoys
6. Realización del docking

---

## Estructura del repositorio

```
PARP1-Virtual-Screening/
├── README.md                   # Este archivo
├── parp1_pipeline.yml          # Entorno conda reproducible
├── scripts/
│   ├── 01_descargar_PARP1_modelo_AF.py  # Descarga del modelo completo de AlphaFold (P09874, v6)
│   ├── 01_descargar_pdbs.py        # Descarga por lotes desde RCSB PDB usando tabla_pdb.csv
│   ├── 02_procesar_pdb.py          # Extracción de la cadena A y el ligando de cada estructura PDB
│   └── 03_contactos.py             # Identificación de interacciones entre ligando y proteína
│   └── 04_lectura.py               # Consenso de las interacciones
│   └── 05_farmacoforo.py           # Interacciones y su relevancia a la hora de la actividad
├── structures/                     
│   └── raw/
│   └── filtered/
│   └── contacts/             
└── info_adicional
    ├── esquema.png               # Esquema del pipeline
    └── tabla_pdb.csv             # Tabla de estructuras PDB y sus ligandos junto a otra información
```

> **Nota**: La carpeta `structures/` está excluida del control de versiones (`.gitignore`) por el tamaño de los archivos PDB. Los scripts de descarga permiten regenerar todas las estructuras localmente.
> **Nota** El archivo IDEAS.md será eliminado una vez terminado el trabajo

---

## Estructuras utilizadas

[Estructuras](Info_adicional/Tabla_pdb.csv)


Todas las estructuras corresponden al **dominio catalítico (CAT)** de PARP1 humana, compuesto por los subdominios HD (661-786) y ART (786-1014), a excepción de la molécula completa predicha usando AlphaFold (PARP1_AF).

Las estructuras correspondientes a los pdb 9ETQ y 9ETR no están disponibles para descargar en formato pdb.

---

## Residuos clave del bolsillo catalítico (REVISAR CON NUEVO PROGRAMA plip)

[Residuos clave](results/farmacoforo.csv)

---

## Instalación y uso

### 1. Clonar el repositorio

```bash
git clone https://github.com/AndoniM12/PARP1-Virtual-Screening.git
cd PARP1-Virtual-Screening
```

### 2. Crear el entorno conda

```bash
conda env create -f parp1_pipeline.yml -n parp1_pipeline
conda activate parp1_pipeline
```

### 3. Descargar las estructuras

**Modelo completo de AlphaFold (PARP1 humana, UniProt P09874):**

```bash
python scripts/01_descargar_PARP1_modelo_AF.py
```

**Estructuras cristalográficas desde RCSB PDB:**

```bash
python scripts/01_descargar_pdbs.py
```

Los PDB descargados se guardan en `structures/raw/`. El archivo `scripts/codigos_pdb.txt` contiene la lista de códigos a descargar (uno por línea, se admiten comentarios con `#`).

### 4. Extraer la cadena A de cada estructura y el ligando que nos interesa

```bash
python scripts/02_procesar_pdb.py
```

Genera versiones `*_chainA.pdb` dentro de `structures/raw/`, listas para su uso en análisis de contactos

### 5. Identificar los contactos entre el ligando y la proteína

```bash
python scripts/03_contactos.py
```

Genera un nuevo directorio `structures/contacts/` con subdirectorios para cada pdb analizado. Dentro de estos directorios se encuentra los archivos .xml para cada código pdb

### 6. Lectura de las interacciones a partir de los archivos xml para cada pdb

```bash
python scripts/04_lectura.py
```

Lee los archivos xml de tadas las estructuras pdb y extrae datos de interes de los contactos detectados en el paso anterior. Archivo `interacciones.csv` creado en el nuevo directorio `reuslts/`

### 7. Generación de farmacoforo

```bash
python scripts/05_farmacoforo.py
```

Crea un consenso teniendo en cuenta las veces que se repite la interacción en cada pdb dando una idea de la importancia de cada interacción y residuo implicado. Genera un nuevo archivo `farmacoforo.csv` en el directorio `reuslts/`


---

## Entorno de trabajo

| Herramienta | Versión | Uso |
|-------------|---------|-----|
| Python | 3.11 | Scripts del pipeline |
| Biopython | 1.87 | Lectura y procesado de PDB |
| NumPy | 2.4.6 | Cálculos numéricos |
| PLIP (Protein-Ligand Interaction Profiles) | 3.0.0 | Identificación automatizada de interacción proteína-ligando|
| UCSF Chimera | — | Inspección y validación visual |
| Pharmit | — | Construcción del farmacóforo y cribado virtual |  - ¿LigandScout como alternativa?
| AlphaFold DB | v6 | Modelo de longitud completa de PARP1 |
| RCSB PDB | — | Estructuras cristalográficas |


---

## Estado actual del proyecto

- [x] Descarga automatizada de estructuras (AlphaFold + PDB)
- [x] Extracción de cadena catalítica
- [x] Análisis de contactos proteína-ligando (plip)
- [ ] Identificación de residuos clave del bolsillo
- [ ] Construcción del farmacóforo de consenso (Pharmit)
- [ ] Validación retrospectiva (activos conocidos + decoys DUD-E)
- [ ] Cribado virtual contra ZINC
- [ ] Análisis y filtrado de hits

---

## Referencias


---

## Autor

**Andoni Moreno Lanceta**  
Trabajo de Fin de Máster — Bioinformática