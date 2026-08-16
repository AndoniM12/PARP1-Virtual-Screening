from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina
import bisect
import os

# ==========================================
# CONFIGURACIÓN
# ==========================================

INPUT_SDF = "results/PHARMIT/top1000_rmsd.sdf"

TOP_RMSD = 1000

OUTPUT_TOP = "results/PHARMIT/top1000_rmsd.sdf"

OUTPUT_FINAL_PATH = "results/DOCKING"
OUTPUT_FINAL = f"{OUTPUT_FINAL_PATH}/representantes_docking.sdf"

CLUSTER_CUTOFF = 0.2


# ==========================================
# Seleccionar top RMSD sin cargar todo en RAM
# ==========================================

def select_top_rmsd(filename, number):

    suppl = Chem.SDMolSupplier(filename)

    top = []

    processed = 0
    counter = 0

    for mol in suppl:

        if mol is None:
            continue

        if not mol.HasProp("rmsd"):
            continue

        processed += 1

        rmsd = float(mol.GetProp("rmsd"))

        # El contador evita comparar objetos Mol cuando hay RMSD iguales
        bisect.insort(top, (rmsd, counter, mol))

        if len(top) > number:
            top.pop()

        counter += 1

        if processed % 100000 == 0:
            print(f"Procesadas {processed:,} moléculas...")

    return [mol for _, _, mol in top]

# ==========================================
# Guardar SDF
# ==========================================

def save_sdf(mols, filename):

    writer = Chem.SDWriter(filename)

    for mol in mols:
        writer.write(mol)

    writer.close()


# ==========================================
# Fingerprints
# ==========================================


def fingerprints(mols):

    generator = AllChem.GetMorganGenerator(
        radius=2,
        fpSize=2048
    )

    fps = []

    for mol in mols:
        fp = generator.GetFingerprint(mol)
        fps.append(fp)

    return fps

# ==========================================
# Clustering
# ==========================================

def clustering(fps):

    distances = []

    for i in range(1, len(fps)):

        similarities = DataStructs.BulkTanimotoSimilarity(
            fps[i],
            fps[:i]
        )

        distances.extend(
            [1 - x for x in similarities]
        )

    clusters = Butina.ClusterData(
        distances,
        len(fps),
        CLUSTER_CUTOFF,
        isDistData=True
    )

    return clusters


# ==========================================
# Seleccionar representantes
# ==========================================

def representatives(mols, clusters):

    selected = []

    for cluster in clusters:

        best = None
        best_rmsd = float("inf")

        for idx in cluster:

            mol = mols[idx]

            rmsd = float(mol.GetProp("rmsd"))

            if rmsd < best_rmsd:

                best = mol
                best_rmsd = rmsd

        selected.append(best)

    return selected


# ==========================================
# MAIN
# ==========================================

def main():

    print("\nSeleccionando TOP RMSD...")

    top1000 = select_top_rmsd(
        INPUT_SDF,
        TOP_RMSD
    )

    print(f"\nMoléculas seleccionadas: {len(top1000)}")

    save_sdf(
        top1000,
        OUTPUT_TOP
    )

    print("\nCalculando fingerprints...")

    fps = fingerprints(top1000)

    print("\nRealizando clustering...")

    clusters = clustering(fps)

    print(f"Clusters obtenidos: {len(clusters)}")

    print("\nSeleccionando representantes...")

    final = representatives(
        top1000,
        clusters
    )

    os.makedirs(
        OUTPUT_FINAL_PATH,
        exist_ok=True
    )

    save_sdf(
        final,
        OUTPUT_FINAL
    )

    print(f"\nMoléculas finales para docking: {len(final)}")


if __name__ == "__main__":
    main()