import os
import pandas as pd
import numpy as np

INPUT_FILE = "results/interacciones.csv"
OUTPUT_FILE = "results/farmacoforo.csv"

def classify_interaction(conservation):

    if conservation >=0.80 :
        return "Obligatory"
    
    elif conservation >= 0.60:
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
            x=("x_pos", "mean"),
            y=("y_pos", "mean"),
            z=("z_pos", "mean"),
            Std_x=("x_pos", "std"),
            Std_y=("y_pos", "std"),
            Std_z=("z_pos", "std"),
            Vector_x=("Vector_x", "mean"),
            Vector_y=("Vector_y", "mean"),
            Vector_z=("Vector_z", "mean"),
            Mean_distance=("Distance", "mean"),
        ).reset_index()
    )

    consensus["Radius"] = np.sqrt(
        consensus["Std_x"]**2 +
        consensus["Std_y"]**2 +
        consensus["Std_z"]**2
        )
    consensus.loc[consensus["PDB"] == 1, "Radius"] = 1


    norm = np.sqrt(
        consensus["Vector_x"]**2 +
        consensus["Vector_y"]**2 +
        consensus["Vector_z"]**2
    )

    consensus["Vector_x"] /= norm
    consensus["Vector_y"] /= norm
    consensus["Vector_z"] /= norm


    consensus["Conservation"] = (
        consensus["PDB"] / total_pdb
    ).round(2)
    consensus["Classification"] = (
        consensus["Conservation"].apply(classify_interaction)
    )
    consensus = consensus.sort_values(
        by=["Conservation"],
        ascending=[False]
    )


    consensus[["x", "y", "z", "Std_x", "Std_y", "Std_z", "Vector_x", "Vector_y", "Vector_z", "Mean_distance", "Radius"]] = (
    consensus[["x", "y", "z", "Std_x", "Std_y", "Std_z", "Vector_x", "Vector_y", "Vector_z", "Mean_distance", "Radius"]]
    .round(3)
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