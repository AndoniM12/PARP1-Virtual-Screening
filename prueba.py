import xml.etree.ElementTree as ET

import os


Input = "structures/contacts/"


for file in os.listdir(Input):

    print(f"Procesando archivo: {file}")
    tree = ET.parse(os.path.join(Input, file, f"{file}_chainA_report.xml"))
    root = tree.getroot()

    binding_site = root.find("bindingsite")

    interactions = binding_site.find("interactions")

    for interaction in interactions:
        print(f"{interaction.tag}: {interaction.text}")

    print("--------------------------------------------------------------------------------------------")
