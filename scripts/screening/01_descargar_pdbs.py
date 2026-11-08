"""
Reúne en data/structures/raw todas las estructuras que necesita el pipeline.

Las estructuras marcadas con Descarga_pdb = "Sí" en tabla_pdb.csv se descargan
del RCSB. Las demás NO están disponibles en formato PDB en el RCSB (devuelve
404), así que vienen incluidas en el repositorio dentro de data/raw/ y este
script solo las copia a su sitio.

Es el caso de 9ETQ (saruparib / AZD5305), que pertenece al Training set y por
tanto entra en el farmacóforo de consenso: sin ella el consenso se construye
sobre una estructura menos y todos los resultados posteriores cambian.
"""

import shutil
import urllib.request
import os
import pandas as pd

TABLE = "data/raw/tabla_pdb.csv"

# Carpeta del repositorio con las estructuras que no se pueden descargar.
BUNDLED_FOLDER = "data/raw"

# Directorio donde se reúnen todas las estructuras de partida.
OUTPUT_FOLDER = "data/structures/raw"

df = pd.read_csv(TABLE)

# Descargables del RCSB.
pdb_codes = df.loc[df["Descarga_pdb"] == "Sí", "PDB"].tolist()

# Incluidas en el repositorio: ni descargables ni generadas de otra forma.
#
# Se excluyen:
#   - PARP1_AF, que tiene su propia descarga desde AlphaFold DB.
#   - las de grupo "-" (8VYH y 9ETR), que el paso 02 no procesa y que el
#     pipeline no usa; estan en la tabla solo como referencia bibliografica.
bundled_codes = df.loc[
    (df["Descarga_pdb"] != "Sí")
    & (df["PDB"] != "PARP1_AF")
    & (df["Grupo de estudio"] != "-"),
    "PDB",
].tolist()

def download_pdb(pdb_code, output_folder):
    """
    Descarga un archivo PDB desde el RCSB PDB y lo guarda en la carpeta especificada.
    """
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    output_path = os.path.join(output_folder, f"{pdb_code}.pdb")
    
    if os.path.exists(output_path):
        print(f"El archivo {output_path} ya existe. Saltando descarga.")
        return
    
    print(f"Descargando {pdb_code} desde {url}\n")
    urllib.request.urlretrieve(url, output_path)
    print(f"Guardado en {output_path}")

def copy_bundled_pdb(pdb_code, output_folder):
    """
    Copia una estructura incluida en el repositorio a la carpeta de trabajo.

    Se usa para las que el RCSB no sirve en formato PDB. No se descarga nada:
    el fichero ya esta versionado en data/raw/.
    """
    source = os.path.join(BUNDLED_FOLDER, f"{pdb_code}.pdb")
    output_path = os.path.join(output_folder, f"{pdb_code}.pdb")

    if os.path.exists(output_path):
        print(f"El archivo {output_path} ya existe. Saltando copia.")
        return True

    if not os.path.exists(source):
        print(
            f"FALTA {pdb_code}: no se puede descargar del RCSB y no esta en "
            f"{source}.\n"
            f"  Esta estructura deberia venir con el repositorio. Sin ella el "
            f"pipeline da resultados distintos."
        )
        return False

    print(f"Copiando {pdb_code} desde {source} (no descargable del RCSB)\n")
    shutil.copyfile(source, output_path)
    print(f"Guardado en {output_path}")
    return True


def download_AF_model(output_folder):

    URL = "https://alphafold.ebi.ac.uk/files/AF-P09874-F1-model_v6.pdb"

    print(f"Descargando estructura de PARP1 (UniProt P09874)...")
    print(f"URL: {URL}")

    OUTPUT_FILENAME = "PARP1_AF_model.pdb"

    urllib.request.urlretrieve(URL, f"{output_folder}/{OUTPUT_FILENAME}")

if __name__ == "__main__":
    # Crear el directorio de salida si no existe
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    # Leer los códigos PDB desde el archivo
    print(f"Se han encontrado {len(pdb_codes)} códigos PDB para descargar.")
    print(f"Y {len(bundled_codes)} incluidos en el repositorio: "
          f"{', '.join(bundled_codes)}\n")

    for code in pdb_codes:
        download_pdb(code, OUTPUT_FOLDER)
        print("------------------------------------------------")

    faltan = []

    for code in bundled_codes:
        if not copy_bundled_pdb(code, OUTPUT_FOLDER):
            faltan.append(code)
        print("------------------------------------------------")

    download_AF_model(OUTPUT_FOLDER)

    print()
    print("=" * 60)
    total = len(pdb_codes) + len(bundled_codes) - len(faltan)
    print(f"Estructuras reunidas en {OUTPUT_FOLDER}: {total}, mas el modelo "
          f"de AlphaFold.")

    if faltan:
        print(f"AVISO  No se han podido reunir: {', '.join(faltan)}")
        print("       Los resultados NO coincidiran con los publicados.")