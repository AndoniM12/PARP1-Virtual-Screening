import pandas as pd

INPUT_FILE = "results/farmacoforo.csv"
OUTPUT_FILE = "results/farmacoforo.pml"

def main():

    df = pd.read_csv(INPUT_FILE)
    
    with open(OUTPUT_FILE, "w") as pml:
        
        pml.write("hide everything\n")
        pml.write("show cartoon\n")
        pml.write("bg_color white\n")

        for n, row in df.iterrows():

            print(row["Feature"])

            feature_name = f"feature_{n}"

            pml.write(
            f"pseudoatom {feature_name}, "
            f"pos=[{row['x']},{row['y']},{row['z']}]\n"
            )

            pml.write(f"show spheres, {feature_name}\n")

            if pd.notna(row["Radius"]):

                pml.write(
                    f"alter {feature_name}, vdw={row['Radius']}\n"
                )

            pml.write("rebuild\n")

if __name__ == "__main__":
    main()