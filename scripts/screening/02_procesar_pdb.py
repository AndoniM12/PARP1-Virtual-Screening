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

INPUT_FOLDER = "data/structures/raw"
OUTPUT_FOLDER = "data/structures/processed"
reference_pdb = "data/structures/processed/7AAA_chainA.pdb"

df = pd.read_csv("data/raw/tabla_pdb.csv")
ligandos_por_pdb = dict(zip(df["PDB"], df["Código_ligando"]))
rmsd_values = []
faltantes = []


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
        pdb_files.remove("7AAA.pdb")

    print(f"Se han encontrado {len(pdb_files)} estructuras para procesar.\n")

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    if not os.path.exists("data/structures/ligands"):
            os.mkdir("data/structures/ligands")

    print(f"Procesando 7AAA...")

    clean_pdb(f"{INPUT_FOLDER}/7AAA.pdb", reference_pdb)

    for _, row in df.iterrows():

        code = row["PDB"]
        ligand = row["Código_ligando"]
        group = row["Grupo de estudio"]

        if group != "-" and code != "PARP1_AF":

            input_path = os.path.join(INPUT_FOLDER, f"{code}.pdb")
            output_path = os.path.join(OUTPUT_FOLDER, f"{code}_chainA.pdb")

            if os.path.exists(output_path):
                print(f"  -> {code} ya procesado. Se omite.")
                continue

            if not os.path.exists(input_path):
                # No se aborta: se avisa y se sigue, para no perder el
                # trabajo ya hecho ni ocultar las demas estructuras que
                # pudieran faltar. El resumen final las lista todas.
                print(f"  -> FALTA {input_path}. Se omite {code}.")
                faltantes.append(code)
                continue

            print(f"Procesando {code}...")

            clean_pdb(input_path, output_path)

            rmsd = align_pdb(reference_pdb, output_path)

            rmsd_values.append(rmsd)

            print(f"  -> Guardado en: {output_path}")

    if faltantes:
        print("\n===== Estructuras que faltan =====")
        print(f"No se han procesado {len(faltantes)}: {', '.join(faltantes)}")
        print("Ejecuta antes scripts/screening/01_descargar_pdbs.py.")
        print("Si alguna pertenece al Training set, el farmacoforo de "
              "consenso saldra distinto y los resultados NO coincidiran "
              "con los publicados.")

    print("\n===== Resumen RMSD =====")
    print(f"Estructuras alineadas contra la referencia : {len(rmsd_values)}")
    print(f"RMSD medio            : {np.mean(rmsd_values):.3f} Å")
    print(f"Desviación estándar   : {np.std(rmsd_values):.3f} Å")
    print(f"RMSD mínimo           : {np.min(rmsd_values):.3f} Å")
    print(f"RMSD máximo           : {np.max(rmsd_values):.3f} Å")