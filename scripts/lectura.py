import os
import xml.etree.ElementTree as ET

import pandas as pd

INPUT_FOLDER = "structures/contacts"
OUTPUT_FILE = "results/interacciones.csv"


def hydrophobic_interactions(binding_site, pdb_code, ligand_name):

    interactions = []

    interactions_node = binding_site.find("interactions")

    if interactions_node is None:
        return interactions
    
    hydrophobic = interactions_node.find("hydrophobic_interactions")

    for interaction in hydrophobic:
     
        interaction_data = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Hydrophobic interaction",
            "Residue": interaction.findtext("restype"),
            "Residue_number": int(interaction.findtext("resnr")),
            "Distance": float(interaction.findtext("dist"))
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
     
        interaction = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Hydrogen bond",
            "Residue": bond.findtext("restype"),
            "Residue_number": int(bond.findtext("resnr")),
            "Distance": float(bond.findtext("dist_d-a"))
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
     
        interaction = {
            "PDB": pdb_code,
            "Ligand": ligand_name,
            "Interaction": "Pi stacking",
            "Residue": stacking.findtext("restype"),
            "Residue_number": int(stacking.findtext("resnr")),
            "Distance": float(stacking.findtext("centdist"))
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

    print(df.head())

if __name__ == "__main__":
    main()