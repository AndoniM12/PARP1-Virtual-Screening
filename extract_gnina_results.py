import os
import csv


INPUT_FOLDER = "results/DOCKING/results_mol2"
OUTPUT_FILE = "gnina_summary.csv"


def parse_gnina_log(file_path):

    ligand = os.path.basename(file_path).replace(".log", "").replace("docked_", "")

    poses = []

    with open(file_path, "r") as f:
        lines = f.readlines()


    start = False

    for line in lines:

        if "mode |  affinity" in line:
            start = True
            continue


        if start:

            values = line.split()

            # buscamos líneas tipo:
            # 1 -7.05 0.64 0.1051 5.140

            if len(values) == 5:

                try:

                    mode = int(values[0])

                    affinity = float(values[1])

                    intramol = float(values[2])

                    cnnscore = float(values[3])

                    cnnaffinity = float(values[4])

                    if affinity <= 0:

                        poses.append({

                            "Ligando": ligand,

                            "Pose": mode,

                            "Vina_score": affinity,

                            "Intramol": intramol,

                            "CNNscore": cnnscore,

                            "CNNaffinity": cnnaffinity

                        })

                        break

                except ValueError:
                    pass


    return poses



all_results = []


for file in sorted(os.listdir(INPUT_FOLDER)):

    if file.endswith(".log"):

        path = os.path.join(
            INPUT_FOLDER,
            file
        )


        results = parse_gnina_log(path)

        all_results.extend(results)



with open(
    OUTPUT_FILE,
    "w",
    newline=""
) as csvfile:


    fieldnames = [

        "Ligando",
        "Pose",
        "Vina_score",
        "Intramol",
        "CNNscore",
        "CNNaffinity"

    ]


    writer = csv.DictWriter(
        csvfile,
        fieldnames=fieldnames
    )


    writer.writeheader()

    writer.writerows(all_results)



print(
    "Archivo generado:",
    OUTPUT_FILE
)


print(
    "Número total de poses:",
    len(all_results)
)