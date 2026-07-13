import pandas as pd

df = pd.read_csv("info_adicional/tabla_pdb.csv")
groups = dict(zip(df["PDB"], df["Grupo de estudio"]))

print(groups)