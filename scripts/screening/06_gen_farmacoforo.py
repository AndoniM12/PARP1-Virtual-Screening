"""
Genera los tres modelos de farmacoforo SIN REFINAR a partir de la tabla de
consenso, en dos formatos:

- .pml  escena de PyMOL para mirar el bolsillo y las esferas
- .json sesion de Pharmit

Estos son los modelos previos al refinamiento. Los que se usan para cribar son
los refinados, que escribe 09_refinamiento_farmacoforo.py en results/pharmit/.

Aqui SI es correcto filtrar por Classification: sobre la tabla de consenso cada
caracteristica tiene una sola clase y no se ha eliminado ni fusionado nada
todavia. Lo que no se puede filtrar asi es la tabla de union que produce el
refinamiento (ver docs/provenance.md).
"""

import json

import pandas as pd

INPUT_FILE = "results/pharmacophore/farmacoforo.csv"
OUTPUT_FILE = "results/pharmacophore/farmacoforo"
REFERENCE_PDB = "data/structures/processed/7AAA_chainA.pdb"
levels = {
    "all" : ["Obligatory", "Recommended", "Optional"],
    "recommended" : ["Obligatory", "Recommended"],
    "core" : ["Obligatory"]
}

# Nombres de caracteristica que entiende Pharmit.
PHARMIT_FEATURE_NAMES = {
    "HBD": "HydrogenDonor",
    "HBA": "HydrogenAcceptor",
    "Hydrophobic": "Hydrophobic",
    "Aromatic": "Aromatic",
}


def write_pharmit_json(filtered_df, path):
    """Escribe la sesion de Pharmit de un modelo sin refinar."""

    points = []

    for _, row in filtered_df.iterrows():

        points.append({
            "name": PHARMIT_FEATURE_NAMES[row["Feature"]],
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
                "z": row["Vector_z"],
            },
            "minsize": "",
            "maxsize": "",
            "selected": True,
        })

    with open(path, "w") as handle:
        json.dump({"points": points}, handle, indent=4)

    return len(points)

def feature_color(feature):

    colors = {
        "HBD": "blue",
        "HBA": "red",
        "Hydrophobic": "yellow",
        "Aromatic": "forest"
    }

    return colors.get(feature, "white")

def feature_transparency(classification):

    transparency = {
        "Obligatory": 0.15,
        "Recommended": 0.35,
        "Optional": 0.65
    }

    return transparency.get(classification, 0.50)

def main():

    df = pd.read_csv(INPUT_FILE)

    for level, categories in levels.items():

        filtered_df = df[df["Classification"].isin(categories)]

        n = write_pharmit_json(filtered_df, f"{OUTPUT_FILE}_{level}.json")
        print(f"{level:<12s} {n:>2d} caracteristicas  ->  "
              f"{OUTPUT_FILE}_{level}.json y _{level}.pml")

        with open(f"{OUTPUT_FILE}_{level}.pml", "w") as pml:

            # Cargar proteína        
            pml.write(f"load {REFERENCE_PDB}, protein\n")
            pml.write("hide everything\n")
            pml.write("show cartoon, protein\n")
            pml.write("color gray90, protein\n")
            pml.write("set cartoon_transparency, 0.65\n")
            pml.write("bg_color white\n")

            # Mostrar bolsillo excluyendo la categoria óptimos
            df_optimal = df[df["Classification"] != "Optional"]
            
            residues = "+".join(
                map(str, sorted(df_optimal["Residue_number"].unique()))
            )

            pml.write(f"select pocket, protein and resi {residues}\n")
            pml.write("show sticks, pocket\n")
            pml.write("color gray70, pocket and elem C\n")
            pml.write("color blue, pocket and elem N\n")
            pml.write("color red, pocket and elem O\n")
            pml.write("color yellow, pocket and elem S\n")
            pml.write("label pocket and name CA, resn + resi\n")
            pml.write("set label_size, 18\n"
                "set label_color, black\n"
                "set label_font_id, 7\n"
                "set label_position, (2,1,1)\n"
                )

            # Dibujar features
            for n, row in filtered_df.iterrows():

                feature_name = f"feature_{n}"

                color = feature_color(row["Feature"])
                transparency = feature_transparency(
                    row["Classification"]
                )

                pml.write(
                    f"pseudoatom {feature_name}, "
                    f"pos=[{row['x']},{row['y']},{row['z']}]\n"
                )
                pml.write(f"hide nonbonded, {feature_name}\n")
                pml.write(f"show spheres, {feature_name}\n")
                pml.write(f"color {color}, {feature_name}\n")

                if pd.notna(row["Radius"]):

                    pml.write(
                        f"alter {feature_name}, vdw={row['Radius']}\n"
                    )

                pml.write(
                    f"set sphere_transparency, "
                    f"{transparency}, {feature_name}\n"
                )

            pml.write("rebuild\n")
            pml.write(    
                "set_view (\\\n"
                "0.295335442,    0.495653808,    0.816764057,\\\n"
                "0.644749045,    0.527470946,   -0.553234816,\\\n"
                "-0.705034792,    0.689999223,   -0.163792834,\\\n"
                "0.000000000,    0.000000000,  -80.127052307,\\\n"
                "13.514476776,   40.709243774,    9.326225281,\\\n"
                "-8592.375000000, 8752.620117188,  -20.000000000 )\n"
            )

            pml.write("set internal_gui, 0")

if __name__ == "__main__":
    main()