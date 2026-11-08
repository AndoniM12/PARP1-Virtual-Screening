from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np


# ============================================================
# ARCHIVOS
# ============================================================

sdf_file = "data/md/ligand_453/ligando_453.sdf"
pdb_file = "results/docking/candidatos_complejo/docked_ligando_453_complex.pdb"


# ============================================================
# CARGAR SDF
# ============================================================

mol = Chem.SDMolSupplier(
    sdf_file,
    removeHs=False
)[0]

if mol is None:
    raise RuntimeError("No se pudo leer el SDF.")

print(f"SDF: {mol.GetNumAtoms()} átomos")


# ============================================================
# EXTRAER ÁTOMOS PESADOS DEL SDF
# ============================================================

sdf_heavy = [
    atom.GetIdx()
    for atom in mol.GetAtoms()
    if atom.GetAtomicNum() > 1
]

print(f"SDF: {len(sdf_heavy)} átomos pesados")


# ============================================================
# LEER LIGANDO DEL PDB
# ============================================================

pdb_heavy = []

with open(pdb_file) as f:
    for line in f:

        if not (
            line.startswith("ATOM")
            or line.startswith("HETATM")
        ):
            continue

        residue = line[17:20].strip()

        if residue != "UNL":
            continue

        element = line[76:78].strip()

        # Si el elemento no está indicado correctamente,
        # usamos la columna correspondiente del PDB.
        if not element:
            element = line[12:14].strip()[0]

        if element.upper() == "H":
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        pdb_heavy.append(
            (element.upper(), np.array([x, y, z]))
        )


print(f"PDB: {len(pdb_heavy)} átomos pesados")


# ============================================================
# COMPROBAR NÚMERO DE ÁTOMOS
# ============================================================

if len(sdf_heavy) != len(pdb_heavy):
    raise RuntimeError(
        f"El número de átomos pesados no coincide: "
        f"SDF={len(sdf_heavy)}, PDB={len(pdb_heavy)}"
    )


# ============================================================
# OBTENER COORDENADAS DEL SDF
# ============================================================

conf = mol.GetConformer()

sdf_coords = []
sdf_elements = []

for idx in sdf_heavy:

    atom = mol.GetAtomWithIdx(idx)

    pos = conf.GetAtomPosition(idx)

    sdf_elements.append(
        atom.GetSymbol().upper()
    )

    sdf_coords.append(
        np.array([pos.x, pos.y, pos.z])
    )


# ============================================================
# COMPROBAR ELEMENTOS Y ORDEN
# ============================================================

print("\nComprobación de elementos:")

for i, (sdf_el, (pdb_el, _)) in enumerate(
    zip(sdf_elements, pdb_heavy),
    start=1
):

    if sdf_el != pdb_el:

        print(
            f"⚠️ Átomo {i}: "
            f"SDF={sdf_el}, PDB={pdb_el}"
        )


# ============================================================
# CALCULAR RMSD DIRECTO
# ============================================================

sdf_coords = np.array(sdf_coords)

pdb_coords = np.array([
    coord
    for _, coord in pdb_heavy
])


differences = sdf_coords - pdb_coords

distances = np.linalg.norm(
    differences,
    axis=1
)

rmsd = np.sqrt(
    np.mean(distances ** 2)
)

max_difference = np.max(distances)

mean_difference = np.mean(distances)


# ============================================================
# RESULTADOS
# ============================================================

print("\n============================================")
print("COMPARACIÓN SDF vs DOCKING PDB")
print("============================================")

print(
    f"RMSD directo:       {rmsd:.6f} Å"
)

print(
    f"Diferencia media:   {mean_difference:.6f} Å"
)

print(
    f"Diferencia máxima:  {max_difference:.6f} Å"
)


# ============================================================
# INTERPRETACIÓN
# ============================================================

if rmsd < 0.01:

    print(
        "\n✓ Las coordenadas son prácticamente idénticas."
    )

elif rmsd < 0.10:

    print(
        "\n✓ La diferencia es muy pequeña y compatible "
        "con redondeo/coordenadas PDB."
    )

elif rmsd < 0.50:

    print(
        "\n⚠️ Existe una diferencia apreciable. "
        "No continuar todavía con la construcción del complejo."
    )

else:

    print(
        "\n✗ Las estructuras no presentan la misma pose. "
        "No utilizar el SDF directamente para construir "
        "las coordenadas del complejo."
    )