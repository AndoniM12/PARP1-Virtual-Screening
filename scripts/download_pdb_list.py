"""
Reads a list of PDB codes from a text file (one code per line)
and downloads each structure directly from the RCSB PDB website,
saving them into the 'structures' folder.
"""

import urllib.request
import os

# Path to the text file containing the PDB codes (one per line)
CODES_FILE = "scripts/pdb_codes.txt"

# Folder where the downloaded files will be saved
OUTPUT_FOLDER = "structures/raw"

def read_pdb_codes(file):
    """
    Reads the text file and returns a list of PDB codes,
    ignoring empty lines, comment lines (starting with #),
    and extra spaces.
    """
    with open(file, "r") as f:
        lines = f.readlines()

    codes = []
    for line in lines:
        line = line.strip()  # remove spaces and newline characters (\n)
        if line and not line.startswith("#"):
            codes.append(line.upper())

    return codes


def download_pdb(pdb_code, output_folder):
    """
    Downloads a single PDB file from RCSB given its code.
    """
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    output_path = os.path.join(output_folder, f"{pdb_code}.pdb")

    print(f"Downloading {pdb_code} from {url}")
    urllib.request.urlretrieve(url, output_path)
    print(f"Saved to {output_path}")


if __name__ == "__main__":
    # Create the output folder if it doesn't exist yet
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Folder created: {OUTPUT_FOLDER}")

    # Read the list of codes and download each one
    pdb_codes = read_pdb_codes(CODES_FILE)
    print(f"Found {len(pdb_codes)} PDB codes to download.\n")

    for code in pdb_codes:
        download_pdb(code, OUTPUT_FOLDER)