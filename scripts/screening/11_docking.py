import os
import subprocess
from rdkit import Chem

# =============================================================================
# Rutas de trabajo
# =============================================================================

# Directorio principal del pipeline de docking
WORK_DIR = "results/DOCKING"

# Biblioteca final de representantes seleccionados tras el clustering
DOCKING_RPR = f"{WORK_DIR}/representantes_docking.sdf"

# Receptor en formato PDBQT
RECEPTOR = f"{WORK_DIR}/7AAA.mol2"

# Directorio donde se almacenarán los ligandos individuales
LIGANDS_SDF = f"{WORK_DIR}/ligands_sdf"

# Directorio donde se almacenarán las poses obtenidas por GNINA
RESULTS_MOL2 = f"{WORK_DIR}/results_mol2"


# =============================================================================
# Definición de la caja de búsqueda
# =============================================================================

# Centro del bolsillo catalítico (Å)
CENTER_X = 8.59
CENTER_Y = 39.35
CENTER_Z = 9.53

# Dimensiones de la caja (Å)
SIZE_X = 22
SIZE_Y = 16
SIZE_Z = 24


# =============================================================================
# Funciones
# =============================================================================

def extract_mols(sdf_input, output_dir):
    """
    Separa un archivo SDF con múltiples moléculas en archivos SDF
    individuales para realizar el docking de cada compuesto por separado.
    """

    os.makedirs(output_dir, exist_ok=True)

    print(f"Extrayendo moléculas de {sdf_input}")

    suppl = Chem.SDMolSupplier(sdf_input)

    count = 0

    for i, molecule in enumerate(suppl):

        if molecule is None:
            continue

        mol_name = f"ligando_{i+1:03d}"

        writer = Chem.SDWriter(
            os.path.join(output_dir, f"{mol_name}.sdf")
        )

        writer.write(molecule)
        writer.close()

        count += 1

    print(f"Moléculas extraídas: {count}")
    print("-------------------------------")


def import_receptor():
    """
    Comprueba si el receptor ya existe en formato PDBQT.
    En caso contrario, lo prepara automáticamente mediante Meeko.
    """

    if os.path.exists(RECEPTOR):
        return

    print("Preparando receptor...")

    receptor_pdb = os.path.join(WORK_DIR, "7AAA_chainA.pdb")

    subprocess.run([
        "cp",
        "structures/processed/7AAA_chainA.pdb",
        receptor_pdb
    ], check=True)

    subprocess.run([
        "mk_prepare_receptor.py",
        "-i", receptor_pdb,
        "-o", RECEPTOR
    ], check=True)

    os.remove(receptor_pdb)

    print("Receptor preparado correctamente mediante Meeko.")
    print("-------------------------------")


def docking(input_file, output_file):
    """
    Ejecuta GNINA para un ligando.

    Para cada compuesto se generan:
        - Un archivo MOL2 con las poses obtenidas.
        - Un archivo LOG con la salida completa de GNINA.
    """

    print(f"Docking: {os.path.basename(input_file)}")

    log_file = output_file.replace(".mol2", ".log")

    with open(log_file, "w") as log:

        subprocess.run([
            "gnina",
            "-r", RECEPTOR,
            "-l", input_file,
            "-o", output_file,
            "--center_x", str(CENTER_X),
            "--center_y", str(CENTER_Y),
            "--center_z", str(CENTER_Z),
            "--size_x", str(SIZE_X),
            "--size_y", str(SIZE_Y),
            "--size_z", str(SIZE_Z),
            "--num_modes", "10",
            "--exhaustiveness", "8",
            "--seed", "1"
        ],
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True)

    return log_file


def extract_docking_results(log_file, ligand_name, csv_file):
    """
    Extrae las puntuaciones de GNINA para todas las poses generadas
    y las almacena en un único archivo CSV.
    """

    start = False

    with open(log_file, "r") as log:

        for line in log:

            line = line.strip()

            if line.startswith("mode"):
                start = True
                continue

            if not start:
                continue

            if line.startswith("-----") or line == "":
                continue

            parts = line.split()

            if len(parts) < 5:
                continue

            if not parts[0].isdigit():
                continue

            pose = int(parts[0])
            affinity = float(parts[1])
            cnn_score = float(parts[3])
            cnn_affinity = float(parts[4])

            with open(csv_file, "a") as csv:

                csv.write(
                    f"{ligand_name},{pose},{affinity},{cnn_score},{cnn_affinity}\n"
                )


# =============================================================================
# Programa principal
# =============================================================================

def main():

    # Preparación del receptor
    import_receptor()

    # Separación de la biblioteca de representantes
    extract_mols(
        DOCKING_RPR,
        LIGANDS_SDF
    )

    # Crear directorio de resultados
    os.makedirs(RESULTS_MOL2, exist_ok=True)

    # Archivo resumen con todas las puntuaciones obtenidas
    results_csv = os.path.join(
        WORK_DIR,
        "results.csv"
    )

    # Crear el CSV e incluir la cabecera
    with open(results_csv, "w") as csv:

        csv.write(
            "Ligando,Pose,Affinity,CNNscore,CNNaffinity\n"
        )

    # Ejecutar el docking para cada ligando
    ligands = sorted(
        os.listdir(LIGANDS_SDF),
        key=lambda f: int(
            f.replace("ligando_", "").replace(".sdf", "")
        )
    )

    for ligand_file in ligands:

        input_file = os.path.join(
            LIGANDS_SDF,
            ligand_file
        )

        output_file = os.path.join(
            RESULTS_MOL2,
            f"docked_{ligand_file[:-4]}.mol2"
        )

        # Ejecutar GNINA
        log_file = docking(
            input_file,
            output_file
        )

        # Extraer las puntuaciones del archivo LOG
        extract_docking_results(
            log_file,
            ligand_file[:-4],
            results_csv
        )

    print("\n==========================================")
    print("Docking molecular finalizado correctamente")
    print(f"Resultados guardados en:\n{results_csv}")
    print("==========================================")


if __name__ == "__main__":
    main()