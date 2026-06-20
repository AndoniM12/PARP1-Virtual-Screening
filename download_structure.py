"""
Downloads the AlphaFold-predicted structure model for human PARP1
(UniProt ID: P09874) from the AlphaFold Protein Structure Database
(https://alphafold.ebi.ac.uk).
"""

import urllib.request
import os

# --- Configuration ---

# UniProt identifier for human PARP1
UNIPROT_ID = "P09874"

# Model version in AlphaFold DB (manually checked on the website (https://alphafold.ebi.ac.uk/entry/P09874) -> currently v6)
VERSION = 6

# Download URL, built following AlphaFold DB's fixed naming pattern
URL = f"https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v{VERSION}.pdb"

# Folder where the downloaded file will be saved
OUTPUT_FOLDER = "structures"

# Name of the output file
OUTPUT_FILENAME = f"AF-{UNIPROT_ID}-F1-model_v{VERSION}.pdb"

# Full path where the file will be saved
OUTPUT_PATH = os.path.join(OUTPUT_FOLDER, OUTPUT_FILENAME)


def download_structure():
    """
    Downloads the PARP1 PDB file from AlphaFold DB and saves it
    in the 'structures' folder. Creates the folder if it doesn't exist.
    """

    # Create the output folder if it doesn't exist yet
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Folder created: {OUTPUT_FOLDER}")

    print(f"Downloading PARP1 structure (UniProt {UNIPROT_ID})...")
    print(f"URL: {URL}")

    # Download the file
    urllib.request.urlretrieve(URL, OUTPUT_PATH)

    print(f"Download complete.\nFile saved at: {OUTPUT_PATH}")


# Entry point of the script: only runs if this file is executed directly
# (not if it gets imported as a module from another script)
if __name__ == "__main__":
    download_structure()