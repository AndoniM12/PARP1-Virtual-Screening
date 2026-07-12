"""
Extrae el contenido de la cadena A de un archivo PDB además de solo la proteína
y el ligando deseado y lo guarda en un nuevo archivo en la carpeta structures/filtered.
"""
import os

import pandas as pd

INPUT_FOLDER = "structures/raw"
OUTPUT_FOLDER = "structures/filtered"


df = pd.read_csv("info_adicional/Tabla_pdb.csv")
ligandos_por_pdb = dict(zip(df["PDB"], df["Código_ligando"]))


def clean_pdb(pdb, output_path):

    with open(pdb, "r") as file:
        for line in file:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22] == "A":
                if line.startswith("ATOM") or (
                    line.startswith("HETATM") and line[17:20].strip() in ligandos_por_pdb.values()
                        ):
                    with open(output_path, "a") as output_file:
                        output_file.write(line)
                        


if __name__ == "__main__":

    pdb_files = os.listdir(INPUT_FOLDER)
    
    if "PARP1_AF_model.pdb" in pdb_files:
        pdb_files.remove("PARP1_AF_model.pdb")
    
    print(f"Se han encontrado {len(pdb_files)} estructuras para procesar.\n")

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    for code in pdb_files:

        input_path = os.path.join(INPUT_FOLDER, f"{code[0:4]}.pdb")
        output_path = os.path.join(OUTPUT_FOLDER, f"{code[0:4]}_chainA.pdb")

        if os.path.exists(output_path):
            print(f"  -> {code} ya procesado. Se omite.")
            continue

        print(f"Procesando {code}...")
        clean_pdb(input_path, output_path)
        print(f"  -> Guardado: {output_path}")