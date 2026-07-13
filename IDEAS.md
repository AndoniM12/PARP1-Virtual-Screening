# TFM - Desarrollo y validación de un pipeline computacional para la identificación de inhibidores de PARP1 a partir de datos públicos

## Objetivo general

Este TFM no pretende descubrir un nuevo inhibidor de PARP1.

El objetivo es desarrollar un pipeline reproducible para cribado virtual y docking molecular que pueda reutilizarse en otras dianas farmacológicas, utilizando PARP1 como caso de estudio debido a la gran cantidad de información estructural y experimental disponible.

---

# Estructura 

1. Introducción

1.1 Cáncer y daño en el ADN
    • Origen del cáncer
    • Tipos de daño en el ADN
    • Importancia de la estabilidad genómica

1.2 Mecanismos de reparación del ADN
    • Reparación de SSB
    • Reparación de DSB
    • Principales vías (BER, NER, MMR, HR y NHEJ)

1.3 La familia PARP
    • Función general
    • Reacción de PARilación
    • Papel del NAD+
    • Diferencias entre los miembros de la familia

1.4 PARP1
    • Organización estructural
    • Dominios funcionales
    • Mecanismo de activación
    • PARilación y automodificación
    • Regulación de la reparación del ADN

1.5 PARP1 como diana terapéutica
    • Relación entre PARP1 y cáncer
    • Concepto de letalidad sintética
    • Mutaciones BRCA1/BRCA2
    • Inhibidores aprobados
    • Mecanismo de acción de los inhibidores
    • Resistencia a los inhibidores
    • Necesidad de desarrollar nuevos compuestos

1.6 Descubrimiento de fármacos asistido por ordenador
    • Drug Discovery in silico
    • Importancia de las estructuras cristalográficas
    • Farmacóforos
    • Virtual Screening
    • Docking molecular

2. Objetivos
   2.1 Objetivo general
   2.2 Objetivos específicos

3. Materiales y métodos
   3.1 Selección de estructuras cristalográficas
   3.2 Obtención de inhibidores conocidos
   3.3 Caracterización del bolsillo catalítico
   3.4 Generación del farmacóforo
   3.5 Validación del farmacóforo
   3.6 Virtual Screening
   3.7 Docking molecular
   3.8 Automatización del pipeline

4. Resultados
   4.1 Selección de la estructura
   4.2 Análisis del bolsillo
   4.3 Farmacóforo obtenido
   4.4 Resultados del cribado virtual
   4.5 Resultados del docking
   4.6 Validación del pipeline

5. Discusión

6. Conclusiones

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
- os
- xml.etree.ElementTree

## Bioinformática estructural

- PyMOL
- PLIP
- AutoDock Vina (o GNINA)

---

# Dudas abiertas

- ¿Validar docking mediante redocking? /  ¿Cómo evaluar el farmacóforo?

Tres grupos de estudio, generación del farmacoforo usando Training set y redocking con Validation set para ver si el farmacóforo creado recoge dichos fármacos

- ¿Qué resolución mínima aceptar?
- ¿Qué scoring utilizar?
- ¿Qué conjunto de compuestos usar para el screening?

ZINC20 / PUBCHEM ...

---

# Próximas tareas

- [ ] Revisar si existe algún tipo de farmacóforo ya modelado
- [x] Crear tabla comparativa de PDB.
- [ ] Descargar ligandos desde ChEMBL.
- [ ] Diseñar la validación del pipeline.

---

# Ideas futuras

- Automatizar completamente el pipeline.

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


## Preguntar por la licencia de LigandScout