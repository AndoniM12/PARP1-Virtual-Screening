"""
Descarga los PDB de la lista de códigos en el archivo codigos_pdb.txt 
y los guarda en la carpeta structures/raw.
"""

import urllib.request
import os
import pandas as pd

# Ruta del archivo de texto que contiene los códigos PDB
df = pd.read_csv("info_adicional/Tabla_pdb.csv")
pdb_codes_dict = dict(zip(df["PDB"],df["Descarga_pdb"]))
pdb_codes = []
for pdb, descarga in pdb_codes.items():
    print(pdb, descarga)
    if descarga == "Sí":
        pdb_codes.append(pdb)