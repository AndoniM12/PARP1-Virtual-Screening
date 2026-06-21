"""
Extracts only chain A from every PDB file in structure/raw
and saves each one as a new file with the suffix '_chainA.pdb'.
"""

from Bio.PDB import PDBParser, PDBIO, Select
import os

INPUT_FOLDER = "structures/raw"
OUTPUT_FOLDER = "structures/chain_A"

class ChainASelector(Select):
    """
    Tells Biopython's PDBIO which atoms to keep when writing
    the new file: only those belonging to chain A.
    """

    def accept_chain(self, chain):
        return chain.id == "A"


def extract_chain_a(pdb_path, output_path):
    """
    Reads a PDB file, keeps only chain A, and writes it to a new file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure/raw", pdb_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, select=ChainASelector())


if __name__ == "__main__":

    pdb_files = os.listdir(INPUT_FOLDER)
    
    print(f"Found {len(pdb_files)} structures to process.\n")

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)
        print(f"Folder created: {OUTPUT_FOLDER}")

    for code in pdb_files:
        input_path = os.path.join(INPUT_FOLDER, f"{code[0:4]}.pdb")
        output_path = os.path.join(OUTPUT_FOLDER, f"{code[0:4]}_chainA.pdb")

        print(f"Processing {code}...")
        extract_chain_a(input_path, output_path)
        print(f"  -> Saved: {output_path}")