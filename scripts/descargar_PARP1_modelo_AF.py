"""
Descarga la estructura predicha por AlphaFold para la proteína humana PARP1
(UniProt ID: P09874) desde la base de datos de estructuras de proteínas AlphaFold
(https://alphafold.ebi.ac.uk).
"""

import urllib.request
import os

# --- Configuration ---

# Identificador UniProt para PARP1 humano
UNIPROT_ID = "P09874"

# URL para descargar el archivo PDB de AlphaFold para PARP1
URL = f"https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v6.pdb"

# Directorio donde se guardará el archivo descargado
OUTPUT_FOLDER = "structures/raw"

# Nombre del archivo de salida (PDB) para la estructura de PARP1
OUTPUT_FILENAME = f"PARP1_AF_model.pdb"

# Ruta completa del archivo de salida
OUTPUT_PATH = os.path.join(OUTPUT_FOLDER, OUTPUT_FILENAME)


def download_structure():
    """
    Descarga el archivo PDB de PARP1 desde la base de datos de AlphaFold y lo guarda
    en la carpeta 'structures'. Crea la carpeta si no existe.
    """

    # Crear la carpeta de salida si no existe
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    print(f"Descargando estructura de PARP1 (UniProt {UNIPROT_ID})...")
    print(f"URL: {URL}")

    # Descargar el archivo PDB y guardarlo en la ruta especificada
    urllib.request.urlretrieve(URL, OUTPUT_PATH)

    print(f"Descarga completada.\nArchivo guardado en: {OUTPUT_PATH}")


# Ejecutar la función de descarga si el script se ejecuta directamente
if __name__ == "__main__":
    download_structure()