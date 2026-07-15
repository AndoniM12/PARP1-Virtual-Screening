import os
import pandas as pd
import numpy as np

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
            ["Interaction", "Residue", "Residue_number", "Feature"]
        ).agg(
            PDB=("PDB", "nunique"),
            Mean_distance=("Distance", "mean"),
            x=("x_pos", "mean"),
            y=("y_pos", "mean"),
            z=("z_pos", "mean"),
            Std_x=("x_pos", "std"),
            Std_y=("y_pos", "std"),
            Std_z=("z_pos", "std")
        ).reset_index()
    )

    consensus["Conservation"] = (
        consensus["PDB"] / total_pdb
    ).round(2)

    consensus["Classification"] = (
        consensus["Conservation"].apply(classify_interaction)
    )

    consensus["Radius"] = np.sqrt(
        consensus["Std_x"]**2 +
        consensus["Std_y"]**2 +
        consensus["Std_z"]**2
    ).round(3)

    consensus = consensus.sort_values(
        by=["Conservation", "Mean_distance"],
        ascending=[False, True]
    )

    consensus[["x", "y", "z", "Mean_distance","Std_x", "Std_y", "Std_z", "Radius"]] = (
    consensus[["x", "y", "z", "Mean_distance","Std_x", "Std_y", "Std_z", "Radius"]]
    .round(2)
    )
    return consensus

def main():

    df = pd.read_csv(INPUT_FILE)

    consensus = consensus_classification(df)

    consensus.to_csv(
        OUTPUT_FILE,
        index = False
    )

    print(consensus)

if __name__ == "__main__":
    main()