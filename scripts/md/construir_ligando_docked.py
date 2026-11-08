from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import re


SDF_FILE = "data/md/ligand_453/ligando_453.sdf"
MOL2_FILE = "data/md/ligand_453/ligando_453.mol2"

print("=" * 60)
print("MAPEO SDF ↔ MOL2 DE GNINA")
print("=" * 60)


# ============================================================
# 1. Cargar SDF
# ============================================================

print("\n1. Cargando SDF...")

supplier = Chem.SDMolSupplier(
    SDF_FILE,
    removeHs=False,
    sanitize=True
)

sdf = supplier[0]

if sdf is None:
    raise RuntimeError("No se pudo cargar el SDF.")

print(f"SDF: {sdf.GetNumAtoms()} átomos")


# ============================================================
# 2. Extraer átomos pesados del SDF
# ============================================================

sdf_heavy_indices = [
    atom.GetIdx()
    for atom in sdf.GetAtoms()
    if atom.GetAtomicNum() != 1
]

print(f"SDF átomos pesados: {len(sdf_heavy_indices)}")


# ============================================================
# 3. Leer MOL2 manualmente
# ============================================================

print("\n2. Cargando MOL2 de GNINA...")

mol2_atoms = []
mol2_bonds = []

section = None

with open(MOL2_FILE) as f:

    for line in f:

        line = line.rstrip()

        if line == "@<TRIPOS>ATOM":
            section = "ATOM"
            continue

        if line == "@<TRIPOS>BOND":
            section = "BOND"
            continue

        if line.startswith("@<TRIPOS>"):
            section = None
            continue

        if section == "ATOM" and line:

            parts = line.split()

            idx = int(parts[0])
            name = parts[1]
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            atom_type = parts[5]

            element = re.match(r"[A-Za-z]+", name).group(0)

            element = element.capitalize()

            mol2_atoms.append({
                "index": idx,
                "name": name,
                "element": element,
                "coords": np.array([x, y, z])
            })

        elif section == "BOND" and line:

            parts = line.split()

            if len(parts) >= 3:

                a = int(parts[1])
                b = int(parts[2])

                mol2_bonds.append((a, b))


print(f"MOL2 átomos: {len(mol2_atoms)}")
print(f"MOL2 enlaces: {len(mol2_bonds)}")


# ============================================================
# 4. Crear molécula RDKit del MOL2 ignorando órdenes de enlace
# ============================================================

print("\n3. Construyendo grafo químico del MOL2...")

mol2_heavy = [
    atom for atom in mol2_atoms
    if atom["element"] != "H"
]

heavy_ids = {
    atom["index"]: i
    for i, atom in enumerate(mol2_heavy)
}

rw = Chem.RWMol()

for atom in mol2_heavy:

    rd_atom = Chem.Atom(atom["element"])

    rw.AddAtom(rd_atom)


for a, b in mol2_bonds:

    # Ignorar enlaces que involucren H
    if a not in heavy_ids or b not in heavy_ids:
        continue

    ia = heavy_ids[a]
    ib = heavy_ids[b]

    rw.AddBond(
        ia,
        ib,
        Chem.BondType.SINGLE
    )


mol2_graph = rw.GetMol()


# ============================================================
# 5. Crear versión del SDF con todos los enlaces como SINGLE
# ============================================================

print("Normalizando conectividad del SDF...")

rw_sdf = Chem.RWMol()

for idx in sdf_heavy_indices:

    atom = sdf.GetAtomWithIdx(idx)

    rd_atom = Chem.Atom(atom.GetSymbol())

    rw_sdf.AddAtom(rd_atom)


sdf_heavy_map = {
    old_idx: new_idx
    for new_idx, old_idx in enumerate(sdf_heavy_indices)
}


for bond in sdf.GetBonds():

    a = bond.GetBeginAtomIdx()
    b = bond.GetEndAtomIdx()

    if a not in sdf_heavy_map or b not in sdf_heavy_map:
        continue

    rw_sdf.AddBond(
        sdf_heavy_map[a],
        sdf_heavy_map[b],
        Chem.BondType.SINGLE
    )


sdf_graph = rw_sdf.GetMol()


# ============================================================
# 6. Buscar correspondencia por ELEMENTO + CONECTIVIDAD
# ============================================================

print("\n4. Buscando correspondencia química...")
print("   Se ignoran los órdenes de enlace y la aromaticidad.")

mcs = rdFMCS.FindMCS(
    [sdf_graph, mol2_graph],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareAny,
    ringMatchesRingOnly=False,
    completeRingsOnly=False,
    matchValences=False,
    timeout=60
)

print(f"\nMCS encontrado: {mcs.numAtoms} átomos")


if mcs.numAtoms != len(sdf_heavy_indices):

    raise RuntimeError(
        f"La correspondencia no es completa: "
        f"{mcs.numAtoms}/{len(sdf_heavy_indices)} átomos."
    )


# ============================================================
# 7. Obtener correspondencia
# ============================================================

mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

sdf_match = sdf_graph.GetSubstructMatch(mcs_mol)
mol2_match = mol2_graph.GetSubstructMatch(mcs_mol)


if len(sdf_match) != len(sdf_heavy_indices):
    raise RuntimeError("No se pudo mapear completamente el SDF.")


if len(mol2_match) != len(mol2_heavy):
    raise RuntimeError("No se pudo mapear completamente el MOL2.")


# ============================================================
# 8. Construir mapa SDF → MOL2
# ============================================================

mapping = {}

for sdf_graph_idx, mol2_graph_idx in zip(
    sdf_match,
    mol2_match
):

    sdf_original_idx = sdf_heavy_indices[sdf_graph_idx]

    mol2_atom = mol2_heavy[mol2_graph_idx]

    mapping[sdf_original_idx] = mol2_atom["index"]


# ============================================================
# 9. Mostrar correspondencia
# ============================================================

print("\n5. Correspondencia SDF → MOL2")
print("-" * 60)

for sdf_idx in sdf_heavy_indices:

    mol2_idx = mapping[sdf_idx]

    sdf_atom = sdf.GetAtomWithIdx(sdf_idx)

    mol2_atom = mol2_atoms[mol2_idx - 1]

    print(
        f"SDF {sdf_idx + 1:2d} "
        f"{sdf_atom.GetSymbol():2s}"
        f"  →  "
        f"MOL2 {mol2_idx:2d} "
        f"{mol2_atom['element']:2s}"
    )


# ============================================================
# 10. Comprobar elementos
# ============================================================

print("\n6. Verificando elementos...")

for sdf_idx, mol2_idx in mapping.items():

    sdf_element = sdf.GetAtomWithIdx(sdf_idx).GetSymbol()
    mol2_element = mol2_atoms[mol2_idx - 1]["element"]

    if sdf_element != mol2_element:

        raise RuntimeError(
            f"Error de elemento: "
            f"SDF {sdf_idx + 1} = {sdf_element}, "
            f"MOL2 {mol2_idx} = {mol2_element}"
        )


print("✓ Todos los elementos coinciden.")


# ============================================================
# 11. Calcular RMSD usando la correspondencia
# ============================================================

print("\n7. Comparando coordenadas...")

sdf_conf = sdf.GetConformer()

sdf_coords = []
mol2_coords = []

for sdf_idx, mol2_idx in mapping.items():

    pos = sdf_conf.GetAtomPosition(sdf_idx)

    sdf_coords.append([
        pos.x,
        pos.y,
        pos.z
    ])

    mol2_coords.append(
        mol2_atoms[mol2_idx - 1]["coords"]
    )


sdf_coords = np.array(sdf_coords)
mol2_coords = np.array(mol2_coords)


# RMSD sin superponer.
# Esto permite comprobar que son dos conformaciones distintas
# y NO modifica las coordenadas del docking.

rmsd = np.sqrt(
    np.mean(
        np.sum(
            (sdf_coords - mol2_coords) ** 2,
            axis=1
        )
    )
)


print(f"RMSD SDF ↔ MOL2 (sin alineamiento): {rmsd:.3f} Å")


print("\n" + "=" * 60)
print("✓ MAPEO COMPLETO")
print("=" * 60)

print(
    "\nLas coordenadas del MOL2 de GNINA pueden utilizarse "
    "como coordenadas definitivas del ligando."
)