"""
Genera el Validation Set para validar el farmacóforo.

El script:

- Lee los ligandos del Validation Set desde tabla_pdb.csv.
- Descarga el ligando ideal desde el Chemical Component Dictionary del PDB.
- Genera múltiples conformaciones mediante RDKit.
- Optimiza cada conformación con MMFF94.
- Guarda todas las conformaciones en un único fichero SDF.

Salida:
    info_adicional/validation_set_db.sdf
"""

import os
from urllib.request import urlopen
from urllib.error import HTTPError

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem


TABLE = "info_adicional/tabla_pdb.csv"
OUTPUT = "info_adicional/validation_set_db.sdf"
NCONF = 100
total = 0


def download_ligand(ligand_code):

    url = f"https://files.rcsb.org/ligands/download/{ligand_code}_ideal.sdf"

    try:
        with urlopen(url) as response:
            return response.read().decode("utf-8")

    except HTTPError:
        print(f"No se pudo descargar {ligand_code}")
        return None

def generate_conformers(sdf_text, pdb_code, writer):

    mol = Chem.MolFromMolBlock(
        sdf_text,
        sanitize=True,
        removeHs=False
    )

    if mol is None:
        print(f"No se pudo leer el ligando {pdb_code}")
        return

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()

    params.randomSeed = 42
    params.pruneRmsThresh = 0.35
    params.useRandomCoords = True

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=NCONF,
        params=params
    )

    if len(conf_ids) == 0:

        print(f"No se pudieron generar conformaciones para {pdb_code}")
        return

    try:
        AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            numThreads=0
        )
    except:
        pass

    mol.SetProp("_Name", pdb_code)

    for conf in conf_ids:

        writer.write(
            mol,
            confId=conf
        )

    print(f"{pdb}: {len(conf_ids)} conformaciones")
    total = len(conf_ids)
    return(total)

if __name__ == "__main__":

    df = pd.read_csv(TABLE)

    validation = df[df["Grupo de estudio"] == "Validation set"]

    if os.path.exists(OUTPUT):
        os.remove(OUTPUT)

    writer = Chem.SDWriter(OUTPUT)

    print(f"Se encontraron {len(validation)} ligandos.\n")

    for _, row in validation.iterrows():

        pdb = row["PDB"]
        ligand = row["Código_ligando"]

        print(f"Procesando {pdb} ({ligand})...")

        sdf = download_ligand(ligand)

        if sdf is None:
            continue

        total += generate_conformers(
            sdf,
            pdb,
            writer
        )

    writer.close()

    print("\n==============================")
    print(f"Ligandos procesados: {len(validation)}")
    print(f"Conformaciones generadas: {total}")
    print(f"Media por ligando: {total/len(validation):.1f}")
    print(f"Archivo generado: {OUTPUT}")
    print("\nValidation Set generado correctamente.")
