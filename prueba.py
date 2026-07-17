import pandas as pd

INPUT_FILE = "results/farmacoforo.csv"
OUTPUT_FILE = "results/farmacoforo"
REFERENCE_PDB = "structures/processed/7AAA_chainA.pdb"
levels = {
    "all" : ["Obligatory", "Recommended", "Optional"],
    "recommended" : ["Obligatory", "Recommended"],
    "core" : ["Obligatory"]
}
df = pd.read_csv(INPUT_FILE)

for level, categories in levels.items():
        
    filtered_df = df[df["Classification"].isin(categories)]

    print(filtered_df)