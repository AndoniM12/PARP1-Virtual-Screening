"""
Prepara las estructuras cristalográficas para el pipeline.

El script:

- selecciona la cadena A,
- conserva únicamente la proteína y el ligando de interés,
- elimina el resto de moléculas,
- alinea todas las estructuras respecto a la referencia 7AAA.pdb
- devuelve información sobre el alineamiento
- extrae el ligando de los Validation set y los guarda en la base de datos sdf

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
ligand_list = "info_adicional/validation_set_db.sdf"

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

def extract_ligand(pdb, ligand):

    with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".sdf",
            delete=False,
        ) as tmp:
            output_sdf = tmp.name

    script = f"""
        load {pdb}, prot

        select lig, resn {ligand}

        save {output_sdf}, lig

        quit
    """

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".pml",
        delete=False
    ) as f:

        f.write(script)
        script_name = f.name

    subprocess.run(
        ["pymol", "-cq", script_name],
        check=True,
        capture_output=True,
        text=True
    )

    with open(output_sdf, "r") as infile, open(ligand_list, "a") as outfile:
        outfile.write(infile.read())

    os.remove(output_sdf)
    os.remove(script_name)

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

    print(f"Procesando 7AAA...")

    clean_pdb(f"{INPUT_FOLDER}/7AAA.pdb", reference_pdb)

    for _, row in df.iterrows():

        code = row["PDB"]
        ligand = row["Código_ligando"]
        group = row["Grupo de estudio"]
        download = row["Descarga_pdb"]

        if download == "Sí":

            input_path = os.path.join(INPUT_FOLDER, f"{code}.pdb")
            output_path = os.path.join(OUTPUT_FOLDER, f"{code}_chainA.pdb")

            if os.path.exists(output_path):
                print(f"  -> {code} ya procesado. Se omite.")
                continue

            print(f"Procesando {code}...")

            clean_pdb(input_path, output_path)

            rmsd = align_pdb(reference_pdb, output_path)

            rmsd_values.append(rmsd)

            if group == "Validation set":

                extract_ligand(output_path,ligand)

            print(f"  -> Guardado en: {output_path}")
    
    print("\n===== Resumen RMSD =====")
    print(f"Estructuras alineadas contra la referencia : {len(rmsd_values)}")
    print(f"RMSD medio            : {np.mean(rmsd_values):.3f} Å")
    print(f"Desviación estándar   : {np.std(rmsd_values):.3f} Å")
    print(f"RMSD mínimo           : {np.min(rmsd_values):.3f} Å")
    print(f"RMSD máximo           : {np.max(rmsd_values):.3f} Å")