from pdbfixer import PDBFixer
from openmm.app import PDBFile

import os

# Entrada: estructura experimental tal como se descargo del PDB.
# Salida: a run/, que esta en .gitignore. La copia de la ejecucion
# original esta versionada en data/md/topology/7AAA_repaired.pdb.
input_pdb = "data/raw/7AAA.pdb"
output_pdb = "run/protein/7AAA_repaired.pdb"

os.makedirs(os.path.dirname(output_pdb), exist_ok=True)

# Cargar estructura experimental
fixer = PDBFixer(filename=input_pdb)

# Obtener las cadenas como lista
chains = list(fixer.topology.chains())

# Mantener únicamente la cadena A
chain_indices_to_remove = [
    i for i, chain in enumerate(chains)
    if chain.id != "A"
]

fixer.removeChains(chain_indices_to_remove)

# No reconstruir residuos completos ausentes.
# GLY660 y SER661 están fuera de la región estructural
# experimental y no deben añadirse.
fixer.findMissingResidues()
fixer.missingResidues = {}

# Identificar átomos pesados faltantes
fixer.findMissingAtoms()

# Reconstruir únicamente los átomos faltantes
fixer.addMissingAtoms()

# Los hidrógenos se añadirán posteriormente mediante pdb2gmx.
PDBFile.writeFile(
    fixer.topology,
    fixer.positions,
    open(output_pdb, "w"),
    keepIds=True
)

print("Estructura reparada correctamente.")
print(f"Archivo generado: {output_pdb}")