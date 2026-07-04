"""
Extrae el contenido de la cadena A de un archivo PDB además de solo la proteína
y el ligando deseado y lo guarda en un nuevo archivo.
"""
import os

INPUT_FOLDER = "structures/raw"
OUTPUT_FOLDER = "structures/chain_A"


ligandos_por_pdb = {

    "7AAA": "",
    "7AAB": "UHB",
    "7AAC": "78P",
    "7AAD": "09L",
    "7KK3": "2YQ",
    "7KK4": "09L",
    "7KK5": "3JD",
    "7KK6": "78P",
    "6BHV": "DQV",
    "9DMC": "DQV",
    "9DMC": "APR"

}


def clean_pdb(pdb, output_path):

    with open(pdb, "r") as file:
        for line in file:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22] == "A":
                if line.startswith("ATOM") or (
                    line.startswith("HETATM") and line[17:20].strip() in ligandos_por_pdb.values()
                        ):
                    with open(output_path, "a") as output_file:
                        output_file.write(line)
                        


if __name__ == "__main__":

    pdb_files = os.listdir(INPUT_FOLDER)
    if "PARP1_AF_model.pdb" in pdb_files:
        pdb_files.remove("PARP1_AF_model.pdb")
    
    print(f"Se han encontrado {len(pdb_files)} estructuras para procesar.\n")

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    for code in pdb_files:
        input_path = os.path.join(INPUT_FOLDER, f"{code[0:4]}.pdb")
        output_path = os.path.join(OUTPUT_FOLDER, f"{code[0:4]}_chainA.pdb")

        print(f"Procesando {code}...")
        clean_pdb(input_path, output_path)
        print(f"  -> Guardado: {output_path}")