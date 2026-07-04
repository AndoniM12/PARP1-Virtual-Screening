# TFM - Desarrollo y validación de un pipeline computacional para la identificación de inhibidores de PARP1 a partir de datos públicos

## Objetivo general

Este TFM no pretende descubrir un nuevo inhibidor de PARP1.

El objetivo es desarrollar un pipeline reproducible para cribado virtual y docking molecular que pueda reutilizarse en otras dianas farmacológicas, utilizando PARP1 como caso de estudio debido a la gran cantidad de información estructural y experimental disponible.

---

# Motivación

¿Por qué PARP1?

- Gran cantidad de estructuras cristalográficas disponibles.
- Numerosos ligandos con actividad experimental.
- Diana farmacológica ampliamente estudiada y de gran interes biológico.
- Excelente candidata para validar un pipeline computacional.

---

# Pregunta de investigación

¿Es posible desarrollar un pipeline reproducible capaz de recuperar y priorizar inhibidores conocidos de PARP1 utilizando únicamente información estructural y bases de datos públicas?

---

# Hipótesis

Un pipeline basado en análisis estructural, generación de farmacóforos, cribado virtual y docking molecular será capaz de recuperar inhibidores conocidos y priorizar nuevos candidatos compatibles con el sitio catalítico de PARP1.

---

# Pipeline propuesto

![Esquema pipeline](Info_adicional/esquema.png)

---

# Decisiones tomadas

## Diana

✔ PARP1

Motivos:

- muchas estructuras
- muchos artículos
- muchos inhibidores
- posibilidad de validación

## Estructuras

[Estructuras](Info_adicional/Tabla_pdb.csv)

---

## Datos

Proteína

- PDB
- AlphaFold

Ligandos

- ChEMBL
- PDB Ligand Expo

Compuestos

- ZINC
- ChEMBL

---

# Herramientas previstas

## Linux

Ubuntu

## Lenguaje

Python

## Librerías

- RDKit
- pandas
- numpy
- matplotlib

## Bioinformática estructural

- PyMOL
- PLIP
- Open Babel
- AutoDock Vina (o GNINA)
- Meeko

---

# Dudas abiertas

- ¿Validar docking mediante redocking?

Validación con las estructuras ya conocidas de los inhibidores ya cristalizados

- ¿Qué resolución mínima aceptar?
- ¿Qué scoring utilizar?
- ¿Cómo evaluar el farmacóforo?
- ¿Qué conjunto de compuestos usar para el screening?

ZINC20 / PUBCHEM ...

---

# Próximas tareas

- [ ] Revisar si existe algún tipo de farmacóforo ya modelado
- [ ] Crear tabla comparativa de PDB.
- [ ] Descargar ligandos desde ChEMBL.
- [ ] Diseñar la validación del pipeline.

---

# Ideas futuras

- Automatizar completamente el pipeline.
- Añadir contenedor Docker.

## Pruebas de herramientas

### OpenPharmacophore

Fecha:
01/07/2026

Instalación:

pip install openpharmacophore

Resultado:
No disponible en PyPI.

Acción:
Instalación desde repositorio GitHub.  (pip install git+https://github.com/UnixJunkie/OpenPharmacophore.git)

Resultado:
No es muy buena idea ya que tiene un monton de dependencias que se tiene que instalar aparte, mejor usar otro programa

## Valorar el uso de RDKit

RDKit podría utilizarse como la base quimioinformática del pipeline, permitiendo trabajar de forma reproducible con las moléculas obtenidas de bases de datos públicas como ChEMBL. Con esta herramienta podrías convertir formatos químicos (SMILES, SDF, MOL), preparar y limpiar los ligandos, calcular propiedades fisicoquímicas (peso molecular, LogP, reglas de Lipinski, donadores y aceptores de puentes de hidrógeno), generar descriptores moleculares y fingerprints para comparar estructuras, agrupar compuestos por similitud, identificar características farmacofóricas (zonas aromáticas, regiones hidrofóbicas, donadores/aceptores de hidrógeno) y ayudar en la selección y filtrado inicial de candidatos antes del docking molecular. Además, al estar integrado con Python, permitiría automatizar estas etapas dentro de un pipeline reproducible y documentado en GitHub.