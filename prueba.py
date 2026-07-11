import xml.etree.ElementTree as ET

import os
import pandas as pd

INPUT_FILE = pd.read_csv("results/interacciones.csv")

total_pdbs = INPUT_FILE["PDB"].nunique()

total_pdb = INPUT_FILE["PDB"].nunique()

consensus = (
    INPUT_FILE.groupby(
        ["Residue", "Residue_number", "Interaction"]
    )["PDB"].nunique().reset_index(name="PDB")
)

consensus["Conservation"] = (
    consensus["PDB"] / total_pdb * 100
)

consensus = consensus.sort_values(
    by=["Conservation"],
    ascending=False
)

print(consensus)
