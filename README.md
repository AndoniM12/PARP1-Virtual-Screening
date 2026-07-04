# PARP1 Virtual Screening Pipeline

Pipeline bioinformático para el cribado virtual de inhibidores de PARP1 (Poli ADP-ribosa polimerasa 1) mediante farmacóforos estructurales. Desarrollado como Trabajo de Fin de Máster (TFM).

---

## Descripción del proyecto

PARP1 es una enzima nuclear implicada en la detección y señalización del daño en el ADN. Su dominio catalítico (CAT), compuesto por los subdominios HD y ART, es la diana terapéutica principal de los inhibidores clínicos de PARP (PARPi), como olaparib, niraparib o veliparib, cuyo mecanismo de acción se basa en la inhibición competitiva frente al cofactor NAD+.

Este proyecto construye un pipeline reproducible y de código abierto para:

1. Descargar y preparar estructuras cristalográficas y modelos predichos de PARP1
2. Analizar las interacciones proteína-ligando en el bolsillo catalítico
3. Construir un farmacóforo de consenso multi-ligando basado en estructura
4. Realizar un cribado virtual contra librerías de moléculas comerciales (ZINC)
5. Validar el farmacóforo mediante validación retrospectiva con activos conocidos y decoys

---

## Estructura del repositorio

```
PARP1-Virtual-Screening/
├── README.md                   # Este archivo
├── project_ideas.md            # Notas y decisiones metodológicas del proyecto
├── parp1_pipeline.yml          # Entorno conda reproducible
├── scripts/
│   ├── PARP1_AF_model_download.py  # Descarga del modelo completo de AlphaFold (P09874, v6)
│   ├── download_pdb_list.py        # Descarga por lotes desde RCSB PDB usando pdb_codes.txt
│   ├── extract_chain_a.py          # Extracción de la cadena A de cada estructura PDB
│   └── pdb_codes.txt               # Lista de códigos PDB utilizados en el análisis
├── structures/
│   └── raw/                        # Estructuras descargadas (no versionadas, ver .gitignore)
│       └── PARP1_AF_model.pdb      # Modelo AlphaFold de PARP1 completa (1014 residuos)
└── images/
    ├── esquema.png                 # Esquema del pipeline
    └── pdb_ligands.png             # Tabla de estructuras PDB y sus ligandos
```

> **Nota**: La carpeta `structures/` está excluida del control de versiones (`.gitignore`) por el tamaño de los archivos PDB. Los scripts de descarga permiten regenerar todas las estructuras localmente.

---

## Estructuras utilizadas

| PDB | Ligando | Tipo | Notas |
|-----|---------|------|-------|
| 7AAA | — | Apo (sin ligando) | Referencia conformacional del CAT |
| 7AAB | EB-47 | Inhibidor experimental | Mismo estudio que 7AAC/7AAD |
| 7AAC | Veliparib | Inhibidor clínico | — |
| 7AAD | Olaparib (09L) | Inhibidor clínico | Referencia principal para el farmacóforo |
| 7KK4 | Olaparib | Inhibidor clínico | Alta calidad, sin huecos en residuos clave |
| 7KK5 | Niraparib | Inhibidor clínico | Alta calidad, sin huecos en residuos clave |
| 6BHV | BAD/DQV | Análogo de NAD+ | HD parcialmente delecionado (intencionado) |
| 9DMC | BAD/DQV + ADP-ribosa | Análogo de NAD+ | His-tag visible en REMARK 465 |

Todas las estructuras corresponden al **dominio catalítico (CAT)** de PARP1 humana, compuesto por los subdominios HD (661-786) y ART (786-1014).

---

## Residuos clave del bolsillo catalítico

Identificados mediante análisis de contactos (`findcontact`) en UCSF Chimera, comparando los inhibidores 09L (olaparib, 7AAD) y DQV (BAD, 6BHV):

| Residuo | Tipo de interacción | Presente en 09L | Presente en DQV |
|---------|---------------------|:-:|:-:|
| His 862 | Puente H | ✓ | ✓ (fuerte, 0.438) |
| Gly 863 | Puente H | ✓ | ✓ |
| Ser 904 | Puente H | ✓ | ✓ |
| Tyr 907 | π-π / aromático | ✓ | ✓ |
| Tyr 896 | π-π / aromático | ✓ | ✓ |
| Tyr 889 | Puente H / aromático | ✓ | ✓ (fuerte) |
| Asp 766 | Puente H | ✓ | — |
| Arg 878 | Puente H | — | ✓ (fuerte, 0.396) |
| Gly 876 | Puente H | — | ✓ |

Los residuos que aparecen en **ambos ligandos** (His862, Gly863, Ser904, Tyr907, Tyr896, Tyr889) constituyen el núcleo del farmacóforo de consenso.

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
python scripts/PARP1_AF_model_download.py
```

**Estructuras cristalográficas desde RCSB PDB:**

```bash
python scripts/download_pdb_list.py
```

Los PDB descargados se guardan en `structures/raw/`. El archivo `scripts/pdb_codes.txt` contiene la lista de códigos a descargar (uno por línea, se admiten comentarios con `#`).

### 4. Extraer la cadena A de cada estructura

```bash
python scripts/extract_chain_a.py
```

Genera versiones `*_chainA.pdb` dentro de `structures/raw/`, listas para su uso en análisis de contactos o construcción del farmacóforo.

---

## Entorno de trabajo

| Herramienta | Versión | Uso |
|-------------|---------|-----|
| Python | 3.11 | Scripts del pipeline |
| Biopython | 1.87 | Lectura y procesado de PDB |
| NumPy | 2.4.6 | Cálculos numéricos |
| UCSF Chimera | — | Visualización y análisis de contactos |
| Pharmit | — | Construcción del farmacóforo y cribado virtual |
| AlphaFold DB | v6 | Modelo de longitud completa de PARP1 |
| RCSB PDB | — | Estructuras cristalográficas |

---

## Estado actual del proyecto

- [x] Descarga automatizada de estructuras (AlphaFold + PDB)
- [x] Extracción de cadena catalítica
- [ ] Análisis de contactos proteína-ligando (Chimera)
- [ ] Identificación de residuos clave del bolsillo
- [ ] Construcción del farmacóforo de consenso (Pharmit)
- [ ] Validación retrospectiva (activos conocidos + decoys DUD-E)
- [ ] Cribado virtual contra ZINC
- [ ] Análisis y filtrado de hits

---

## Referencias

- Conceição et al. (2025). *PARP1: A comprehensive review of its mechanisms, therapeutic implications and emerging cancer treatments*. BBA - Reviews on Cancer, 1880, 189282.
- Koes & Camacho (2014). *Pharmit: interactive exploration of chemical space*. Nucleic Acids Research.
- Suskiewicz et al. (2023). Updated protein domain annotation of the PARP protein family. Nucleic Acids Research.

---

## Autor

**Andoni Moreno Lanceta**  
Trabajo de Fin de Máster — Bioinformática