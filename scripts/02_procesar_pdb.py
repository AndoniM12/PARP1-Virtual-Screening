"""
Prepara las estructuras cristalográficas para el pipeline.

El script:

- selecciona la cadena A,
- conserva únicamente la proteína y el ligando de interés,
- elimina el resto de moléculas,
- alinea todas las estructuras respecto a la referencia 7AAA.pdb
- devuelve información sobre el alineamiento

Las estructuras procesadas se guardan en structures/processed.
"""


import os
import pandas as pd
import subprocess
import tempfile
import numpy as np
import re

INPUT_FOLDER = "structures/raw"
OUTPUT_FOLDER = "structures/processed"
reference_pdb = "structures/processed/7AAA_chainA.pdb"

df = pd.read_csv("info_adicional/tabla_pdb.csv")
ligandos_por_pdb = dict(zip(df["PDB"], df["Código_ligando"]))
rmsd_values = []


def clean_pdb(input_pdb, output_pdb):

    with open(input_pdb) as infile, open(output_pdb, "w") as outfile:

        for line in infile:

            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == "A":

                if line.startswith("ATOM") or (
                    line.startswith("HETATM") and
                    line[17:20].strip() in ligandos_por_pdb.values()
                ):

                    outfile.write(line)
                        
def align_pdb(reference_pdb, pdb):

        script = f"""
    load {reference_pdb}, ref
    load {pdb}, mob

    super mob, ref

    save {pdb}, mob

    quit
    """

        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".pml",
            delete=False
        ) as f:

            f.write(script)
            script_name = f.name

        result = subprocess.run(
            ["pymol", "-cq", script_name],
            check=True,
            capture_output=True,
            text=True
        )

        os.remove(script_name)

        match = re.search(r"Executive:\s*RMSD\s*=\s*(\d+\.\d+)", result.stdout)        
        
        return float(match.group(1))

if __name__ == "__main__":

    pdb_files = os.listdir(INPUT_FOLDER)
    
    if "PARP1_AF_model.pdb" in pdb_files:
        pdb_files.remove("PARP1_AF_model.pdb")

    if "7AAA" in pdb_files:
        pdb_files.remove("PARP1_AF_model.pdb")

    print(f"Se han encontrado {len(pdb_files)} estructuras para procesar.\n")

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    clean_pdb(f"{INPUT_FOLDER}/7AAA.pdb", reference_pdb)

    for code in pdb_files:

        input_path = os.path.join(INPUT_FOLDER, f"{code[:4]}.pdb")
        output_path = os.path.join(OUTPUT_FOLDER, f"{code[:4]}_chainA.pdb")

        if os.path.exists(output_path):
            print(f"  -> {code} ya procesado. Se omite.")
            continue

        print(f"Procesando {code}...")

        clean_pdb(input_path, output_path)

        rmsd = align_pdb(reference_pdb, os.path.join(OUTPUT_FOLDER, f"{code[:4]}_chainA.pdb"))

        rmsd_values.append(rmsd)
        
        print(f"  -> Guardado: {output_path}")
    
    print("\n===== Resumen RMSD =====")
    print(f"Estructuras alineadas contra la referencia : {len(rmsd_values)}")
    print(f"RMSD medio            : {np.mean(rmsd_values):.3f} Å")
    print(f"Desviación estándar   : {np.std(rmsd_values):.3f} Å")
    print(f"RMSD mínimo           : {np.min(rmsd_values):.3f} Å")
    print(f"RMSD máximo           : {np.max(rmsd_values):.3f} Å")