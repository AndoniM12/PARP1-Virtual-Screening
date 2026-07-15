"""
Este script calcula los contactos de las estructuras PDB filtradas utilizando la herramienta PLIP.
Se asume que los archivos PDB filtrados se encuentran en la carpeta "structures/filtered" 
y los resultados se guardarán en la carpeta "structures/contacts" solo manteniendo los archivos XML.
"""

import os
import subprocess 
import pandas as pd  

INPUT_FOLDER = "structures/processed"
OUTPUT_FOLDER = "structures/contacts"

df = pd.read_csv("info_adicional/tabla_pdb.csv")
groups = dict(zip(df["PDB"], df["Grupo de estudio"]))
        
def plip_process (pdb,input_path, output_path):

    print(f"Calculando contactos para {pdb}")

    cmd = [
        "plip",
        "-f", str(input_path),
        "-o", str(output_path),
        "-x",
        "-t"
    ]

    subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )
    
    for archivo in os.listdir(output_path):
                if not archivo.endswith(".xml"):
                    os.remove(os.path.join(output_path, archivo))

def main ():

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    all_pdb = os.listdir(INPUT_FOLDER)
    training_pdb = []

    for file in all_pdb:
        if groups.get(file[:4]) == "Training set":
            training_pdb.append(file)

    print(f"Se han encontrado {len(training_pdb)} estructuras para analizar:\n")

    for code in training_pdb:

        input_path = os.path.join(INPUT_FOLDER, code)
        output_path = os.path.join(OUTPUT_FOLDER, code[:4])

        if os.path.exists(output_path):
            print(f"  -> {code[:4]} ya procesado. Se omite.")
            continue

        os.makedirs(output_path, exist_ok=True)

        plip_process(code[:4], input_path, output_path)

if __name__ == "__main__":
    main()