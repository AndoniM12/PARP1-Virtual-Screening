"""
Descarga los PDB de la lista de códigos en el archivo codigos_pdb.txt 
y los guarda en la carpeta structures/raw.
"""

import urllib.request
import os

# Ruta del archivo de texto que contiene los códigos PDB
CODES_FILE = "scripts/codigos_pdb.txt"

# Directorio donde se guardarán los archivos PDB descargados
OUTPUT_FOLDER = "structures/raw"

def read_pdb_codes(file):
    """
    Lee los códigos PDB desde un archivo de texto y 
    devuelve una lista de códigos.
    """
    with open(file, "r") as f:
        lines = f.readlines()

    codes = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith("#"):
            codes.append(line.upper())

    return codes


def download_pdb(pdb_code, output_folder):
    """
    Descarga un archivo PDB desde el RCSB PDB y lo guarda en la carpeta especificada.
    """
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    output_path = os.path.join(output_folder, f"{pdb_code}.pdb")

    print(f"Descargando {pdb_code} desde {url}")
    urllib.request.urlretrieve(url, output_path)
    print(f"Guardado en {output_path}")


if __name__ == "__main__":
    # Crear el directorio de salida si no existe
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    # Leer los códigos PDB desde el archivo
    pdb_codes = read_pdb_codes(CODES_FILE)
    print(f"Se han encontrado {len(pdb_codes)} códigos PDB para descargar.\n")

    for code in pdb_codes:
        download_pdb(code, OUTPUT_FOLDER)