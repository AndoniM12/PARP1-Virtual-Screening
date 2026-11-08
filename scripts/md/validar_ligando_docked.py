import re
import numpy as np


import os

# Valida lo que se acaba de generar en run/. Si no existe, recurre a
# la copia versionada de la ejecucion original en data/md/.
VERSIONADO = "data/md/ligand_453/gromacs"
GENERADO = "run/ligand/gromacs"


def elegir(nombre):
    ruta = os.path.join(GENERADO, nombre)
    return ruta if os.path.exists(ruta) else os.path.join(VERSIONADO, nombre)


GRO_FILE = elegir("ligand_453_docked.gro")
ITP_FILE = elegir("ligand_453_LIG.itp")

print(f"GRO: {GRO_FILE}")
print(f"ITP: {ITP_FILE}")


# ============================================================
# Leer GRO
# ============================================================

with open(GRO_FILE) as f:
    lines = f.readlines()

n_atoms = int(lines[1].strip())

coords = {}

for line in lines[2:2 + n_atoms]:

    atom_number = int(line[15:20])

    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])

    coords[atom_number] = np.array([x, y, z])


# ============================================================
# Leer enlaces del ITP
# ============================================================

bonds = []

section = None

with open(ITP_FILE) as f:

    for line in f:

        line = line.strip()

        if not line or line.startswith(";"):
            continue

        if line.startswith("["):
            section = line.strip("[]").strip()
            continue

        if section == "bonds":

            parts = line.split()

            if len(parts) >= 2:

                a = int(parts[0])
                b = int(parts[1])

                bonds.append((a, b))


# ============================================================
# Identificar enlaces H-X
# ============================================================

print("=" * 60)
print("VALIDACIÓN GEOMÉTRICA DEL LIGANDO DOCKEADO")
print("=" * 60)

print(f"\nÁtomos: {n_atoms}")
print(f"Enlaces: {len(bonds)}")

print("\nEnlaces que contienen H:")
print("-" * 60)


h_bonds = []

for a, b in bonds:

    # En nuestro GRO los H son 33-53
    if a >= 33 or b >= 33:

        h = a if a >= 33 else b
        heavy = b if a >= 33 else a

        distance = np.linalg.norm(
            coords[h] - coords[heavy]
        )

        h_bonds.append(distance)

        print(
            f"H {h:2d} - átomo {heavy:2d}: "
            f"{distance:.3f} nm "
            f"({distance * 10:.2f} Å)"
        )


# ============================================================
# Estadísticas
# ============================================================

if h_bonds:

    h_bonds = np.array(h_bonds)

    print("\n" + "-" * 60)
    print(
        f"Distancia H-X media: "
        f"{h_bonds.mean():.3f} nm "
        f"({h_bonds.mean() * 10:.2f} Å)"
    )

    print(
        f"Distancia mínima: "
        f"{h_bonds.min():.3f} nm "
        f"({h_bonds.min() * 10:.2f} Å)"
    )

    print(
        f"Distancia máxima: "
        f"{h_bonds.max():.3f} nm "
        f"({h_bonds.max() * 10:.2f} Å)"
    )


# ============================================================
# Criterio de validación
# ============================================================

# Enlaces covalentes X-H suelen estar aproximadamente
# entre 0.09 y 0.11 nm en esta representación.

bad = [
    d for d in h_bonds
    if d < 0.07 or d > 0.13
]


print("\n" + "=" * 60)

if bad:

    print(
        f"⚠️ Se encontraron {len(bad)} "
        f"enlaces H-X fuera del rango esperado."
    )

    print(
        "No debemos montar todavía el complejo."
    )

else:

    print(
        "✓ GEOMETRÍA DE LOS HIDRÓGENOS VALIDADA"
    )

    print(
        "Las distancias H-X son compatibles "
        "con enlaces covalentes."
    )

print("=" * 60)