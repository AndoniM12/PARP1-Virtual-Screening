import os
import pandas as pd

INPUT_FILE = "results/interacciones.csv"
OUTPUT_FILE = "results/farmacoforo.csv"

def classify_interaction(conservation):

    if conservation >= 90:
        return "Obligatory"
    
    elif conservation >= 70:
        return "Recommended"
    else:
        return "Optional"

def consensus_classification(df):

    total_pdb = df["PDB"].nunique()

    consensus = (
        df.groupby(
            ["Residue", "Residue_number", "Interaction"]
        )["PDB"].nunique().reset_index(name="PDB")
    )

    consensus["Conservation"] = (
        consensus["PDB"] / total_pdb * 100
    ).round(2)

    consensus["Classification"] = consensus["Conservation"].apply(classify_interaction)

    consensus = consensus.sort_values(
        by=["Conservation"],
        ascending=False
    )
    
    return consensus

def main():

    df = pd.read_csv(INPUT_FILE)

    consensus = consensus_classification(df)

    consensus.to_csv(
        OUTPUT_FILE
    )

    print(consensus)

if __name__ == "__main__":
    main()