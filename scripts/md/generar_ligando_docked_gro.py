import numpy as np
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import rdFMCS


import os

SDF_FILE = "data/md/ligand_453/ligando_453.sdf"
MOL2_FILE = "data/md/ligand_453/ligando_453.mol2"

# Salida a run/, que esta en .gitignore. La copia de la ejecucion
# original esta versionada en data/md/ligand_453/gromacs/.
GRO_OUTPUT = "run/ligand/gromacs/ligand_453_docked.gro"

os.makedirs(os.path.dirname(GRO_OUTPUT), exist_ok=True)


# ============================================================
# LEER MOL2
# ============================================================

def read_mol2(path):

    atoms = []
    bonds = []
    section = None

    atomic_numbers = {
        "C": 6,
        "N": 7,
        "O": 8,
        "S": 16,
        "H": 1
    }

    with open(path) as f:

        for line in f:

            line = line.strip()

            if line == "@<TRIPOS>ATOM":
                section = "atoms"
                continue

            if line == "@<TRIPOS>BOND":
                section = "bonds"
                continue

            if line.startswith("@<TRIPOS>"):
                section = None
                continue

            if section == "atoms" and line:

                p = line.split()

                # El elemento se obtiene del tipo Tripos
                element = p[1][0].upper()

                atoms.append({
                    "id": int(p[0]),
                    "element": element,
                    "x": float(p[2]),
                    "y": float(p[3]),
                    "z": float(p[4])
                })

            elif section == "bonds" and line:

                p = line.split()

                bonds.append(
                    (
                        int(p[1]) - 1,
                        int(p[2]) - 1
                    )
                )

    return atoms, bonds


# ============================================================
# CONSTRUIR GRAFO RDKit
# ============================================================

def build_graph(elements, bonds):

    atomic_numbers = {
        "C": 6,
        "N": 7,
        "O": 8,
        "S": 16,
        "H": 1
    }

    mol = Chem.RWMol()

    for element in elements:

        mol.AddAtom(
            Chem.Atom(
                atomic_numbers[element]
            )
        )

    for a, b in bonds:

        mol.AddBond(
            a,
            b,
            Chem.BondType.SINGLE
        )

    return mol.GetMol()


# ============================================================
# INICIO
# ============================================================

print("=" * 70)
print("CONSTRUCCIÓN DEL LIGANDO DOCKEADO")
print("=" * 70)


# ============================================================
# OPENFF
# ============================================================

print("\nCargando molécula OpenFF...")

ligand = Molecule.from_file(
    SDF_FILE,
    file_format="SDF"
)

print(f"OpenFF: {ligand.n_atoms} átomos")

if ligand.n_atoms != 53:

    raise RuntimeError(
        "La molécula OpenFF no tiene los 53 átomos esperados."
    )


# ============================================================
# MOL2
# ============================================================

print("\nLeyendo MOL2 de GNINA...")

mol2_atoms, mol2_bonds = read_mol2(
    MOL2_FILE
)

print(f"MOL2: {len(mol2_atoms)} átomos")


# ============================================================
# SEPARAR PESADOS
# ============================================================

heavy_indices = [
    i
    for i, atom in enumerate(ligand.atoms)
    if atom.atomic_number != 1
]

print(
    f"Átomos pesados OpenFF: {len(heavy_indices)}"
)

if len(heavy_indices) != 32:

    raise RuntimeError(
        "No hay exactamente 32 átomos pesados."
    )


# ============================================================
# GRAFO OPENFF PESADO
# ============================================================

elements_openff = [
    ligand.atoms[i].symbol
    for i in heavy_indices
]

heavy_map_old_to_new = {
    old: new
    for new, old in enumerate(heavy_indices)
}


bonds_openff = []

for bond in ligand.bonds:

    a = bond.atom1_index
    b = bond.atom2_index

    if a in heavy_map_old_to_new and b in heavy_map_old_to_new:

        bonds_openff.append(
            (
                heavy_map_old_to_new[a],
                heavy_map_old_to_new[b]
            )
        )


openff_graph = build_graph(
    elements_openff,
    bonds_openff
)


# ============================================================
# GRAFO MOL2
# ============================================================

elements_mol2 = [
    atom["element"]
    for atom in mol2_atoms
]

mol2_graph = build_graph(
    elements_mol2,
    mol2_bonds
)


# ============================================================
# MCS
# ============================================================

print("\nBuscando correspondencia OpenFF ↔ GNINA...")

mcs = rdFMCS.FindMCS(
    [
        openff_graph,
        mol2_graph
    ],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareAny,
    ringMatchesRingOnly=False,
    completeRingsOnly=False,
    timeout=60
)

print(
    f"MCS: {mcs.numAtoms} átomos"
)

if mcs.numAtoms != 32:

    raise RuntimeError(
        "No se encontró correspondencia completa de los átomos pesados."
    )


mcs_mol = Chem.MolFromSmarts(
    mcs.smartsString
)

match_openff = openff_graph.GetSubstructMatch(
    mcs_mol
)

match_mol2 = mol2_graph.GetSubstructMatch(
    mcs_mol
)

mapping = dict(
    zip(
        match_openff,
        match_mol2
    )
)

print("✓ Correspondencia completa.")


# ============================================================
# COORDENADAS GNINA
# ============================================================

mol2_coords = np.array([
    [
        atom["x"],
        atom["y"],
        atom["z"]
    ]
    for atom in mol2_atoms
])


# ============================================================
# COORDENADAS OPENFF
# ============================================================

if ligand.conformers is None:

    raise RuntimeError(
        "OpenFF no contiene conformación."
    )

openff_coords = np.array(
    ligand.conformers[0].magnitude
)


# ============================================================
# POSICIONES DOCKEADAS
# ============================================================

docked_coords = np.zeros(
    (53, 3)
)


# ------------------------------------------------------------
# ÁTOMOS PESADOS
# ------------------------------------------------------------

for openff_graph_idx, mol2_idx in mapping.items():

    original_idx = heavy_indices[
        openff_graph_idx
    ]

    docked_coords[
        original_idx
    ] = mol2_coords[
        mol2_idx
    ]


# ============================================================
# TRANSFORMACIÓN PARA LOS H
# ============================================================

P = []
Q = []

for openff_graph_idx, mol2_idx in mapping.items():

    original_idx = heavy_indices[
        openff_graph_idx
    ]

    P.append(
        openff_coords[
            original_idx
        ]
    )

    Q.append(
        mol2_coords[
            mol2_idx
        ]
    )


P = np.array(P)
Q = np.array(Q)


centroid_P = P.mean(axis=0)
centroid_Q = Q.mean(axis=0)

P0 = P - centroid_P
Q0 = Q - centroid_Q

H = P0.T @ Q0

U, S, Vt = np.linalg.svd(H)

R = Vt.T @ U.T

if np.linalg.det(R) < 0:

    Vt[-1, :] *= -1

    R = Vt.T @ U.T


transformed = (
    (openff_coords - centroid_P)
    @ R.T
    + centroid_Q
)


# ============================================================
# HIDRÓGENOS
# ============================================================

print("\nColocando hidrógenos...")

for i, atom in enumerate(ligand.atoms):

    if atom.atomic_number != 1:
        continue

    neighbors = list(
        atom.bonds
    )

    if len(neighbors) != 1:

        raise RuntimeError(
            f"H {i+1} no tiene exactamente un enlace."
        )

    bond = neighbors[0]

    if bond.atom1_index == i:

        heavy_idx = bond.atom2_index

    else:

        heavy_idx = bond.atom1_index


    # Vector H-X original
    vector = (
        openff_coords[i]
        - openff_coords[heavy_idx]
    )

    norm = np.linalg.norm(vector)

    if norm < 1e-8:

        raise RuntimeError(
            f"Vector H-X inválido para H {i+1}."
        )

    vector /= norm


    # Longitud covalente aproximada
    heavy_element = ligand.atoms[
        heavy_idx
    ].symbol

    if heavy_element == "C":

        distance = 1.09

    elif heavy_element == "N":

        distance = 1.01

    elif heavy_element == "O":

        distance = 0.96

    elif heavy_element == "S":

        distance = 1.34

    else:

        distance = 1.09


    # El vector se transforma con la misma
    # rotación obtenida de los átomos pesados
    vector_docked = vector @ R.T

    vector_docked /= np.linalg.norm(
        vector_docked
    )


    docked_coords[i] = (
        docked_coords[heavy_idx]
        + vector_docked * distance
    )


# ============================================================
# VALIDACIÓN
# ============================================================

print("\nValidando distancias H-X...")

distances = []

for i, atom in enumerate(ligand.atoms):

    if atom.atomic_number != 1:
        continue

    bond = list(atom.bonds)[0]

    if bond.atom1_index == i:

        heavy_idx = bond.atom2_index

    else:

        heavy_idx = bond.atom1_index

    distance = np.linalg.norm(
        docked_coords[i]
        - docked_coords[heavy_idx]
    )

    distances.append(
        distance
    )

    print(
        f"H {i+1:2d} - "
        f"{ligand.atoms[heavy_idx].symbol}"
        f"{heavy_idx+1:2d}: "
        f"{distance:.3f} Å"
    )


distances = np.array(distances)

print("\n" + "-" * 70)

print(
    f"Media: {distances.mean():.3f} Å"
)

print(
    f"Mínima: {distances.min():.3f} Å"
)

print(
    f"Máxima: {distances.max():.3f} Å"
)


if np.any(distances < 0.85) or np.any(distances > 1.45):

    raise RuntimeError(
        "Hay enlaces H-X fuera del rango razonable."
    )


print(
    "\n✓ GEOMETRÍA DE HIDRÓGENOS VALIDADA"
)


# ============================================================
# GENERAR GRO
# ============================================================

with open(GRO_OUTPUT, "w") as f:

    f.write(
        "Ligand 453 - GNINA docking pose\n"
    )

    f.write(
        f"{ligand.n_atoms:5d}\n"
    )

    for i, atom in enumerate(ligand.atoms):

        element = atom.symbol

        x, y, z = (
            docked_coords[i] / 10.0
        )

        f.write(
            f"{1:5d}"
            f"{'LIG':<5s}"
            f"{element:>5s}"
            f"{i+1:5d}"
            f"{x:8.3f}"
            f"{y:8.3f}"
            f"{z:8.3f}\n"
        )

    f.write(
        "  10.00000  10.00000  10.00000\n"
    )


print("\n" + "=" * 70)
print("✓ LIGANDO DOCKEADO GENERADO")
print("=" * 70)

print(
    f"\nArchivo:\n{GRO_OUTPUT}"
)

print(
    f"\nÁtomos: {ligand.n_atoms}"
)