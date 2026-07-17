import os
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np

INPUT_FOLDER = "structures/contacts"
OUTPUT_FILE = "results/interacciones.csv"


def hydrophobic_interactions(binding_site, pdb_code, ligand_name):

    interactions = []

    interactions_node = binding_site.find("interactions")

    if interactions_node is None:
        return interactions
    
    hydrophobic = interactions_node.find("hydrophobic_interactions")


    for interaction in hydrophobic:
    
        ligcoo = interaction.find("ligcoo")
        protcoo = interaction.find("protcoo")
       
        # Vector de interaccion de puente de hidrógeno
        vector_x= float(ligcoo.findtext("x")) - float(protcoo.findtext("x"))
        vector_y= float(ligcoo.findtext("y")) - float(protcoo.findtext("y"))
        vector_z= float(ligcoo.findtext("z")) - float(protcoo.findtext("z"))

        norm = np.sqrt(vector_x**2 + vector_y**2 + vector_z**2)
        vector_x /= norm
        vector_y /= norm
        vector_z /= norm


        interaction_data = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Hydrophobic interaction",
            "Residue": interaction.findtext("restype"),
            "Residue_number": int(interaction.findtext("resnr")),
            "Distance": float(interaction.findtext("dist")),
            "Feature": "Hydrophobic",
            "x_pos": float(ligcoo.findtext("x")),
            "y_pos": float(ligcoo.findtext("y")),
            "z_pos": float(ligcoo.findtext("z")),
            "Vector_x": vector_x.round(3),
            "Vector_y": vector_y.round(3),
            "Vector_z": vector_z.round(3),
        }

        interactions.append(interaction_data)

    return interactions

def hydrogen_bonds(binding_site, pdb_code, ligand_name):

    interactions = []

    interactions_node = binding_site.find("interactions")

    if interactions_node is None:
        return interactions
    
    hbonds = interactions_node.find("hydrogen_bonds")

    for bond in hbonds:
    
        ligcoo = bond.find("ligcoo")
        protcoo = bond.find("protcoo")
        protisdon = bond.findtext("protisdon")
        
        # Vector de interaccion de puente de hidrógeno
        vector_x= float(ligcoo.findtext("x")) - float(protcoo.findtext("x"))
        vector_y= float(ligcoo.findtext("y")) - float(protcoo.findtext("y"))
        vector_z= float(ligcoo.findtext("z")) - float(protcoo.findtext("z"))

        norm = np.sqrt(vector_x**2 + vector_y**2 + vector_z**2)
        vector_x /= norm
        vector_y /= norm
        vector_z /= norm

        
        # Verificación de carácter donador o aceptor
        if protisdon == "true":
            protisdon = "HBA"
        else:
            protisdon = "HBD"

        interaction = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Hydrogen bond",
            "Residue": bond.findtext("restype"),
            "Residue_number": int(bond.findtext("resnr")),
            "Distance": float(bond.findtext("dist_d-a")),
            "Feature": protisdon ,
            "x_pos": float(ligcoo.findtext("x")),
            "y_pos": float(ligcoo.findtext("y")),
            "z_pos": float(ligcoo.findtext("z")),
            "Vector_x": vector_x.round(3),
            "Vector_y": vector_y.round(3),
            "Vector_z": vector_z.round(3)
        }

        interactions.append(interaction)

    return interactions

def pi_stacks(binding_site, pdb_code, ligand_name):

    interactions = []

    interactions_node = binding_site.find("interactions")

    if interactions_node is None:
        return interactions
    
    pi_stacks = interactions_node.find("pi_stacks")

    for stacking in pi_stacks:
        
        ligcoo = stacking.find("ligcoo")
        protcoo = stacking.find("protcoo")
       
        # Vector de interaccion de puente de hidrógeno
        vector_x= float(ligcoo.findtext("x")) - float(protcoo.findtext("x"))
        vector_y= float(ligcoo.findtext("y")) - float(protcoo.findtext("y"))
        vector_z= float(ligcoo.findtext("z")) - float(protcoo.findtext("z"))

        norm = np.sqrt(vector_x**2 + vector_y**2 + vector_z**2)
        vector_x /= norm
        vector_y /= norm
        vector_z /= norm

        interaction = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Pi stacking",
            "Residue": stacking.findtext("restype"),
            "Residue_number": int(stacking.findtext("resnr")),
            "Distance": float(stacking.findtext("centdist")),
            "Feature": "Aromatic",
            "x_pos": float(ligcoo.findtext("x")),
            "y_pos": float(ligcoo.findtext("y")),
            "z_pos": float(ligcoo.findtext("z")),
            "Vector_x": vector_x.round(3),
            "Vector_y": vector_y.round(3),
            "Vector_z": vector_z.round(3),
        }

        interactions.append(interaction)

    return interactions

def parse_report(xml_file):

    tree = ET.parse(xml_file)
    root = tree.getroot()

    binding_site = root.find("bindingsite")

    identifiers = binding_site.find("identifiers")

    pdb_code = root.find("pdbid").text[:4]
    ligand_name = identifiers.findtext("longname")

    interactions = []


    interactions.extend(
        hydrophobic_interactions(
            binding_site,
            pdb_code,
            ligand_name
        )
    )
    
    interactions.extend(
        hydrogen_bonds(
            binding_site,
            pdb_code,
            ligand_name
        )
    )

    interactions.extend(
        pi_stacks(
            binding_site,
            pdb_code,
            ligand_name
        )
    )

    return interactions

def main():

    all_interactions = []

    for folder in sorted(os.listdir(INPUT_FOLDER)):

        xml_file = os.path.join(
            INPUT_FOLDER,
            folder,
            f"{folder}_chainA_report.xml"
        )

        if not os.path.exists(xml_file):
            continue

        print(f"Leyendo {folder}")

        interactions = parse_report(xml_file)

        all_interactions.extend(interactions)

    df = pd.DataFrame(all_interactions)

    if not os.path.exists("results"):

        os.makedirs("results", exist_ok=True)

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Se han guardado {len(df)} interacciones en {OUTPUT_FILE} correspondientes a {df['PDB'].nunique()} estructuras.")

if __name__ == "__main__":
    main()