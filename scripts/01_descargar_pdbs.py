"""
Descarga los PDB de la lista de códigos en el archivo tabla_pdb.csv
y los guarda en la carpeta structures/raw.
"""

import urllib.request
import os
import pandas as pd

# Selección de IDs con posibilidad de descarga a través de PDB
df = pd.read_csv("info_adicional/tabla_pdb.csv")
pdb_codes_dict = dict(zip(df["PDB"],df["Descarga_pdb"]))
pdb_codes = []

for pdb, descarga in pdb_codes_dict.items():

    if descarga == "Sí":
        pdb_codes.append(pdb)

# Directorio donde se guardarán los archivos PDB descargados
OUTPUT_FOLDER = "structures/raw"

def download_pdb(pdb_code, output_folder):
    """
    Descarga un archivo PDB desde el RCSB PDB y lo guarda en la carpeta especificada.
    """
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    output_path = os.path.join(output_folder, f"{pdb_code}.pdb")
    
    if os.path.exists(output_path):
        print(f"El archivo {output_path} ya existe. Saltando descarga.")
        return
    
    print(f"Descargando {pdb_code} desde {url}"\n)
    urllib.request.urlretrieve(url, output_path)
    print(f"Guardado en {output_path}")

if __name__ == "__main__":
    # Crear el directorio de salida si no existe
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    # Leer los códigos PDB desde el archivo
    print(f"Se han encontrado {len(pdb_codes)} códigos PDB para descargar.\n")

    for code in pdb_codes:
        download_pdb(code, OUTPUT_FOLDER)
        print("------------------------------------------------")