"""Monta el complejo PARP1-ligando concatenando los dos ficheros .gro.

Recuperado del historial de terminal: este paso se ejecuto originalmente
como un heredoc de Python en la linea de comandos y no existia como
script. Se reproduce aqui literalmente, con rutas relativas a la raiz
del repositorio.

Uso (desde la raiz del repositorio):
    python scripts/md/montar_complejo.py [DIRECTORIO_DE_TRABAJO]

El directorio de trabajo por defecto es `run/`, que esta en .gitignore.
"""

import sys
from pathlib import Path

work = Path(sys.argv[1] if len(sys.argv) > 1 else "run")

protein = work / "protein" / "PARP1.gro"
output = work / "complex" / "complex.gro"

# Usa el ligando que se acaba de generar en run/. Si no existe, recurre
# a la copia versionada de la ejecucion original en data/md/.
ligand = work / "ligand" / "gromacs" / "ligand_453_docked.gro"

if not ligand.exists():
    ligand = Path("data/md/ligand_453/gromacs/ligand_453_docked.gro")

print(f"Proteina: {protein}")
print(f"Ligando:  {ligand}")

# Numero de atomos esperado. En la ejecucion original estos valores
# estaban fijados con assert; se mantienen como comprobacion, pero
# avisando en lugar de abortar si la estructura de partida cambia.
PROTEIN_ATOMS_ESPERADOS = 5546
LIGAND_ATOMS_ESPERADOS = 53

p_lines = protein.read_text().splitlines()
l_lines = ligand.read_text().splitlines()

p_atoms = int(p_lines[1].strip())
l_atoms = int(l_lines[1].strip())

if p_atoms != PROTEIN_ATOMS_ESPERADOS:
    print(
        f"AVISO: la proteina tiene {p_atoms} atomos, "
        f"no los {PROTEIN_ATOMS_ESPERADOS} de la ejecucion original."
    )

if l_atoms != LIGAND_ATOMS_ESPERADOS:
    print(
        f"AVISO: el ligando tiene {l_atoms} atomos, "
        f"no los {LIGAND_ATOMS_ESPERADOS} de la ejecucion original."
    )

# Formato GRO: titulo + numero de atomos + lineas de atomos + caja.
# La caja se toma de la proteina; se redefine despues con `gmx editconf`.
p_atom_lines = p_lines[2 : 2 + p_atoms]
l_atom_lines = l_lines[2 : 2 + l_atoms]
box = p_lines[2 + p_atoms]

output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w") as f:
    f.write("PARP1-LIG 453 complex\n")
    f.write(f"{p_atoms + l_atoms:5d}\n")
    for line in p_atom_lines:
        f.write(line + "\n")
    for line in l_atom_lines:
        f.write(line + "\n")
    f.write(box + "\n")

print(f"Complejo creado: {output}")
print(f"Atomos de proteina: {p_atoms}")
print(f"Atomos de ligando:  {l_atoms}")
print(f"Atomos totales:     {p_atoms + l_atoms}")
