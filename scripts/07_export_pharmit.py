import os
import pandas as pd
import json

INPUT_FILE = "results/farmacoforo.csv"
OUTPUT_FILE = "results/farmacoforo"
df = pd.read_csv(INPUT_FILE)

levels = {
    "all" : ["Obligatory", "Recommended", "Optional"],
    "recommended" : ["Obligatory", "Recommended"],
    "core" : ["Obligatory"]
}

interactions = {
    "HBD" : "HydrogenDonor",
    "HBA" : "HydrogenAcceptor",
    "Hydrophobic" : "Hydrophobic",
    "Aromatic" : "Aromatic"
}

def export_farmaophore(df):
    
    for level, categories in levels.items():

        filtered_df = df[df["Classification"].isin(categories)]

        session = {
                "points": []
                }
        
        for index, row in filtered_df.iterrows():
            
            point = {
                "name": interactions[row["Feature"]],
                "hasvec": True,
                "x": row["x"],
                "y": row["y"],
                "z": row["z"],
                "radius": row["Radius"],
                "enabled": True,
                "vector_on": 1,
                "svector": {
                    "x": row["Vector_x"],
                    "y": row["Vector_y"],
                    "z": row["Vector_z"]
                },
                "minsize": "",
                "maxsize": "",
                "selected": True
            }

            session["points"].append(point)

        with open(f"{OUTPUT_FILE}_{level}.json", "w") as f:
            json.dump(session, f, indent=4)

if __name__ == "__main__":
    export_farmaophore(df)