"""
09_refine_pharmacophore.py

Refina automáticamente los farmacóforos antes de utilizarlos en Pharmit.

Reglas implementadas:

1. Eliminar regiones hidrofóbicas redundantes con regiones aromáticas.
2. Fusionar regiones hidrofóbicas cercanas.
3. Aplicar radios mínimos.
4. Aplicar radios máximos.

Autor:
Andoni Moreno Lanceta
"""

import json
import math
import os

# ==========================================================
# PARÁMETROS
# ==========================================================

INPUT_FOLDER = "results"
OUTPUT_FOLDER = "results/PHARMIT"

HYDROPHOBIC = "Hydrophobic"
AROMATIC = "Aromatic"

# Solapamiento aromático-hidrofóbico
OVERLAP_FACTOR = 0.90

# Fusión de hidrofóbicas
MERGE_FACTOR = 1.20

MIN_RADIUS = {
    "Hydrophobic": 2.20,
    "Aromatic": 1.75,
    "HydrogenDonor": 0.60,
    "HydrogenAcceptor": 0.60,
}

MAX_RADIUS = {
    "Hydrophobic": 3.00,
    "Aromatic": 2.75,
    "HydrogenDonor": 1.20,
    "HydrogenAcceptor": 1.20,
}

# ==========================================================


def distance(p1, p2):

    return math.sqrt(
        (p1["x"] - p2["x"]) ** 2 +
        (p1["y"] - p2["y"]) ** 2 +
        (p1["z"] - p2["z"]) ** 2
    )


# ==========================================================


def remove_hydrophobic_inside_aromatic(points):

    refined = []

    aromatics = [
        p for p in points
        if p["name"] == AROMATIC
    ]

    removed = 0

    for point in points:

        if point["name"] != HYDROPHOBIC:
            refined.append(point)
            continue

        delete = False

        for aromatic in aromatics:

            d = distance(point, aromatic)

            limit = OVERLAP_FACTOR * (
                aromatic["radius"] +
                point["radius"]
            )

            if d < limit:
                delete = True
                removed += 1
                break

        if not delete:
            refined.append(point)

    return refined, removed


# ==========================================================


def merge_hydrophobics(points):

    hydrophobic = [
        p for p in points
        if p["name"] == HYDROPHOBIC
    ]

    others = [
        p for p in points
        if p["name"] != HYDROPHOBIC
    ]

    merged = []
    visited = set()

    fused = 0

    for i in range(len(hydrophobic)):

        if i in visited:
            continue

        # Construimos un grupo mediante búsqueda (BFS)
        queue = [i]
        group_idx = []

        while queue:

            current_idx = queue.pop(0)

            if current_idx in visited:
                continue

            visited.add(current_idx)
            group_idx.append(current_idx)

            current = hydrophobic[current_idx]

            for j in range(len(hydrophobic)):

                if j in visited:
                    continue

                candidate = hydrophobic[j]

                d = distance(current, candidate)

                limit = MERGE_FACTOR * (
                    current["radius"] +
                    candidate["radius"]
                )

                if d < limit:
                    queue.append(j)

        group = [hydrophobic[k] for k in group_idx]

        fused += len(group) - 1

        if len(group) == 1:

            merged.append(group[0])

        else:

            new_point = group[0].copy()

            # Centroide
            new_point["x"] = sum(p["x"] for p in group) / len(group)
            new_point["y"] = sum(p["y"] for p in group) / len(group)
            new_point["z"] = sum(p["z"] for p in group) / len(group)

            # Radio que cubra toda la región
            new_radius = max(
                distance(new_point, p) + p["radius"]
                for p in group
            )

            new_point["radius"] = round(new_radius, 3)

            merged.append(new_point)

    return others + merged, fused


# ==========================================================


def adjust_radius(points):

    for point in points:

        name = point["name"]

        if name in MIN_RADIUS:

            point["radius"] = max(
                point["radius"],
                MIN_RADIUS[name]
            )

        if name in MAX_RADIUS:

            point["radius"] = min(
                point["radius"],
                MAX_RADIUS[name]
            )

    return points


# ==========================================================


if __name__ == "__main__":

    os.makedirs(
        OUTPUT_FOLDER,
        exist_ok=True
    )

    pharmacophores = sorted([
        f for f in os.listdir(INPUT_FOLDER)
        if f.endswith(".json")
    ])

    print(f"\nSe encontraron {len(pharmacophores)} farmacóforos.\n")

    for file in pharmacophores:

        print("=" * 60)
        print(f"Procesando {file}")

        with open(
            os.path.join(INPUT_FOLDER, file),
            "r"
        ) as f:

            pharmacophore = json.load(f)

        points = pharmacophore["points"]

        original = len(points)

        points, removed = remove_hydrophobic_inside_aromatic(points)

        after_remove = len(points)

        points, fused = merge_hydrophobics(points)

        after_merge = len(points)

        points = adjust_radius(points)

        pharmacophore["points"] = points

        output_name = file.replace(
            ".json",
            "_refined.json"
        )

        with open(
            os.path.join(
                OUTPUT_FOLDER,
                output_name
            ),
            "w"
        ) as f:

            json.dump(
                pharmacophore,
                f,
                indent=4
            )

        print(f"Features iniciales                : {original}")
        print(f"Hidrofóbicas eliminadas           : {removed}")
        print(f"Hidrofóbicas fusionadas           : {fused}")
        print(f"Features finales                  : {after_merge}")

    print("\nRefinamiento completado correctamente.")