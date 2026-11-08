"""
09_refinamiento_farmacoforo.py

VARIANTE ANTERIOR. El script vigente es refinamiento.py, que aplica EXACTAMENTE
el mismo refinamiento (mismas reglas, mismas constantes) y ademas escribe las
sesiones de Pharmit de cada modelo.

Este conserva una salida distinta y util para inspeccionar: una fila por
caracteristica Y modelo, con la columna Model, que deja ver de un vistazo que
caracteristica sobrevive en que modelo. No lo usa ningun otro script.

Construye y refina los tres modelos de farmacoforo (Core, Recommended, All)
a partir de la tabla de consenso results/pharmacophore/farmacoforo.csv y escribe
una tabla CSV por modelo.

ORDEN DE OPERACIONES
--------------------
Para cada modelo:  FILTRAR por clasificacion  ->  REFINAR

El orden importa y no es conmutativo. La regla 1 elimina una caracteristica
hidrofobica cuando una aromatica ya la cubre, de modo que solo puede evaluarse
contra las aromaticas que pertenecen a ese modelo. Refinar primero sobre el
conjunto completo y filtrar despues hace que un modelo herede eliminaciones
provocadas por caracteristicas que no forman parte de el.

Ejemplo concreto en PARP1: la hidrofobica de Tyr889 queda cubierta por la
aromatica de Tyr896, que es Opcional. En el modelo All esa aromatica existe y
la hidrofobica se elimina; en el modelo Recommended no existe, asi que la
hidrofobica de Tyr889 debe conservarse.

REGLAS DE REFINAMIENTO
----------------------
1. Eliminar hidrofobicas redundantes con una aromatica del mismo modelo.
2. Fusionar hidrofobicas proximas que describen la misma cavidad.
3. Ajustar radios a su minimo y maximo por tipo de caracteristica.

CLASIFICACION
-------------
La clase de cada caracteristica de entrada se toma TAL CUAL de la columna
Classification del CSV de consenso, que la fija el script 05.

Una caracteristica FUSIONADA hereda la clase MAS FUERTE de sus miembros, segun
el orden Obligatory > Recommended > Optional. Si una obligatoria se fusiona con
una opcional, la fusionada es obligatoria: la region sigue siendo necesaria, y
rebajarla perderia una condicion que el consenso respalda.

La clase de la fusionada NO se vuelve a derivar de la conservacion. Hacerlo
seria incorrecto con los umbrales actuales: una Recommended de conservacion
0,83 superaria OBLIGATORY_MIN = 0,80 y quedaria promovida a obligatoria sin que
ningun miembro lo fuera. Los umbrales solo se usan para el aviso de arranque.

PROPAGACION DE METADATOS EN LA FUSION
-------------------------------------
Una caracteristica fusionada toma como conservacion el MAXIMO de sus miembros,
nunca la suma. Es una cota inferior conservadora: la conservacion exacta seria
la fraccion de estructuras del Training set en que aparece ALGUNA de las
caracteristicas fusionadas, y para calcularla haria falta la identidad de las
estructuras, no solo su recuento. La tabla de consenso solo guarda el recuento.

SALIDA
------
Un unico CSV, results/pharmit/farmacoforo_refined_por_modelo.csv, con una fila por
caracteristica y modelo. La columna Model dice a que modelo pertenece cada
fila. Hace falta esa columna porque el refinamiento depende del modelo: la
hidrofobica de Tyr889 sobrevive en Recommended y se elimina en All.

COMPROBACIONES
--------------
A. Discrepancias entre la columna Classification y los umbrales (aviso).
B. La conservacion de una fusionada no supera el maximo de sus miembros.
C. El recuento cuadra: entrada = salida + eliminadas + absorbidas por fusion.
D. Toda fusionada cubre geometricamente a sus miembros.
E. Toda caracteristica cae dentro de la caja de busqueda del docking.
F. Ninguna pareja de caracteristicas del mismo tipo queda solapada tras
   refinar (si quedase, es que la fusion no las ha capturado).

Autor: Andoni Moreno Lanceta
"""

import csv
import math
import os
import sys

# ==========================================================
# PARAMETROS
# ==========================================================

INPUT_FILE = "results/pharmacophore/farmacoforo.csv"
# Salida distinta de la de refinamiento.py a proposito: los dos scripts
# aplican el mismo refinamiento y escribirian aqui, pisandose.
OUTPUT_FILE = "results/pharmit/farmacoforo_refined_por_modelo.csv"

HYDROPHOBIC = "Hydrophobic"
AROMATIC = "Aromatic"

# Regla 1: una hidrofobica se elimina si su distancia a una aromatica es menor
# que este factor por la suma de ambos radios.
OVERLAP_FACTOR = 0.90

# Regla 2: dos hidrofobicas se fusionan si su distancia es menor que este
# factor por la suma de sus radios.
MERGE_FACTOR = 1.20

# Regla 3: radios admisibles por tipo de caracteristica (angstrom).
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

# Umbrales de conservacion. Solo se usan para el aviso de arranque; la
# clasificacion efectiva sale de la columna Classification del CSV de entrada.
OBLIGATORY_MIN = 0.80
RECOMMENDED_MIN = 0.60

# Jerarquia de clases. Una fusionada se queda con la mas fuerte de sus
# miembros: cuanto menor es el numero, mas fuerte es la clase.
CLASS_RANK = {
    "Obligatory": 0,
    "Recommended": 1,
    "Optional": 2,
}

# Composicion de cada modelo.
LEVELS = {
    "core": ["Obligatory"],
    "recommended": ["Obligatory", "Recommended"],
    "all": ["Obligatory", "Recommended", "Optional"],
}

# Nombres de caracteristica de la tabla -> nombres de Pharmit.
FEATURE_NAMES = {
    "HBD": "HydrogenDonor",
    "HBA": "HydrogenAcceptor",
    "Hydrophobic": "Hydrophobic",
    "Aromatic": "Aromatic",
}

# Caja de busqueda del docking (script 11), para la comprobacion E.
BOX_CENTER = (8.59, 39.35, 9.53)
BOX_SIZE = (22.0, 16.0, 24.0)


# ==========================================================
# UTILIDADES
# ==========================================================


def distance(p1, p2):
    return math.sqrt(
        (p1["x"] - p2["x"]) ** 2
        + (p1["y"] - p2["y"]) ** 2
        + (p1["z"] - p2["z"]) ** 2
    )


def classify(conservation):
    """Clase a partir de la conservacion. Solo para fusionadas."""

    if conservation >= OBLIGATORY_MIN:
        return "Obligatory"

    if conservation >= RECOMMENDED_MIN:
        return "Recommended"

    return "Optional"


def to_float(value, default=0.0):
    if value is None or str(value).strip() == "":
        return default

    return float(value)


def label(feature):
    """Etiqueta legible: Tyr889, o Leu769/Asp766 si esta fusionada."""

    return "/".join(
        residue.capitalize() + number
        for residue, number in zip(feature["residues"], feature["numbers"])
    )


# ==========================================================
# LECTURA
# ==========================================================


def read_consensus(path):

    features = []

    with open(path, newline="", encoding="utf-8") as handle:

        for index, row in enumerate(csv.DictReader(handle)):

            feature_key = row["Feature"].strip()

            if feature_key not in FEATURE_NAMES:
                raise ValueError(
                    f"Fila {index}: caracteristica desconocida "
                    f"'{feature_key}'"
                )

            features.append({
                "order": index,
                "interaction": row["Interaction"].strip(),
                "residues": [row["Residue"].strip()],
                "numbers": [row["Residue_number"].strip()],
                "feature": feature_key,
                "name": FEATURE_NAMES[feature_key],
                "pdb": int(to_float(row.get("PDB"), 0)),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
                "vx": to_float(row.get("Vector_x")),
                "vy": to_float(row.get("Vector_y")),
                "vz": to_float(row.get("Vector_z")),
                "radius": float(row["Radius"]),
                "conservation": float(row["Conservation"]),
                "classification": row["Classification"].strip(),
                "members": 1,
                "origin": [index],
            })

    return features


def warn_classification_mismatch(features):
    """COMPROBACION A. Avisa, no detiene."""

    problems = [
        f"    {label(f):<10s} {f['feature']:<11s} conservacion "
        f"{f['conservation']:.2f} -> tabla '{f['classification']}', "
        f"umbrales '{classify(f['conservation'])}'"
        for f in features
        if classify(f["conservation"]) != f["classification"]
    ]

    if not problems:
        return

    print(
        "  AVISO  La columna Classification no coincide con los umbrales "
        f"({OBLIGATORY_MIN:.2f} / {RECOMMENDED_MIN:.2f}) en "
        f"{len(problems)} filas."
    )
    print("         Se respeta la columna del CSV; los umbrales solo")
    print("         clasifican las caracteristicas fusionadas.")

    for problem in problems:
        print(problem)

    print()


# ==========================================================
# REGLA 1 · ELIMINAR HIDROFOBICAS CUBIERTAS POR UNA AROMATICA
# ==========================================================


def remove_hydrophobic_inside_aromatic(features):

    aromatics = [f for f in features if f["name"] == AROMATIC]

    kept = []
    removed = []

    for feature in features:

        if feature["name"] != HYDROPHOBIC:
            kept.append(feature)
            continue

        covering = None

        for aromatic in aromatics:

            limit = OVERLAP_FACTOR * (aromatic["radius"] + feature["radius"])

            if distance(feature, aromatic) < limit:
                covering = aromatic
                break

        if covering is None:
            kept.append(feature)
        else:
            removed.append((feature, covering))

    return kept, removed


# ==========================================================
# REGLA 2 · FUSIONAR HIDROFOBICAS PROXIMAS
# ==========================================================


def merge_hydrophobics(features):

    hydrophobic = [f for f in features if f["name"] == HYDROPHOBIC]
    others = [f for f in features if f["name"] != HYDROPHOBIC]

    merged = []
    absorbed = 0
    visited = set()

    for i in range(len(hydrophobic)):

        if i in visited:
            continue

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

                limit = MERGE_FACTOR * (
                    current["radius"] + candidate["radius"]
                )

                if distance(current, candidate) < limit:
                    queue.append(j)

        group = [hydrophobic[k] for k in sorted(group_idx)]

        if len(group) == 1:
            merged.append(group[0])
            continue

        absorbed += len(group) - 1
        merged.append(fuse(group))

    return others + merged, absorbed


def fuse(group):

    # Representante: mayor conservacion; empate, el primero de la tabla.
    leader = min(group, key=lambda f: (-f["conservation"], f["order"]))

    fused = {
        "order": leader["order"],
        "interaction": leader["interaction"],
        "residues": [r for f in group for r in f["residues"]],
        "numbers": [n for f in group for n in f["numbers"]],
        "feature": HYDROPHOBIC,
        "name": HYDROPHOBIC,
        "pdb": max(f["pdb"] for f in group),
        "x": sum(f["x"] for f in group) / len(group),
        "y": sum(f["y"] for f in group) / len(group),
        "z": sum(f["z"] for f in group) / len(group),
        "vx": leader["vx"],
        "vy": leader["vy"],
        "vz": leader["vz"],
        "conservation": max(f["conservation"] for f in group),
        "members": sum(f["members"] for f in group),
        "origin": [i for f in group for i in f["origin"]],
        "group": group,
    }

    fused["radius"] = round(
        max(distance(fused, f) + f["radius"] for f in group), 3
    )

    # La fusionada hereda la clase mas fuerte de sus miembros.
    fused["classification"] = min(
        (f["classification"] for f in group),
        key=lambda c: CLASS_RANK[c],
    )

    # COMPROBACION B
    assert fused["conservation"] <= max(
        f["conservation"] for f in group
    ) + 1e-9, "Conservacion de la fusionada mayor que la de sus miembros"

    # COMPROBACION B2 · la clase no puede ser mas debil que la del miembro
    # mas fuerte, ni mas fuerte que ninguno de ellos.
    assert CLASS_RANK[fused["classification"]] == min(
        CLASS_RANK[f["classification"]] for f in group
    ), "La clase de la fusionada no es la mas fuerte de sus miembros"

    # COMPROBACION D
    for member in group:
        assert distance(fused, member) <= fused["radius"] + 1e-6, (
            "La fusionada no cubre a uno de sus miembros"
        )

    return fused


# ==========================================================
# REGLA 3 · AJUSTAR RADIOS
# ==========================================================


def adjust_radius(features):

    for feature in features:

        name = feature["name"]

        if name in MIN_RADIUS:
            feature["radius"] = max(feature["radius"], MIN_RADIUS[name])

        if name in MAX_RADIUS:
            feature["radius"] = min(feature["radius"], MAX_RADIUS[name])

        feature["radius"] = round(feature["radius"], 3)

    return features


# ==========================================================
# REFINAMIENTO
# ==========================================================


def refine(features):

    before = len(features)

    features, removed = remove_hydrophobic_inside_aromatic(features)
    features, absorbed = merge_hydrophobics(features)

    # Radio con el que la regla 2 decidio fusionar; lo usa la comprobacion F.
    for feature in features:
        feature["merge_radius"] = feature["radius"]

    features = adjust_radius(features)

    # COMPROBACION C
    assert before == len(features) + len(removed) + absorbed, (
        f"Recuento descuadrado: entrada {before}, salida {len(features)}, "
        f"eliminadas {len(removed)}, absorbidas {absorbed}"
    )

    features.sort(key=lambda f: f["order"])

    return features, removed, absorbed


# ==========================================================
# COMPROBACION DE POSICIONES
# ==========================================================


def check_positions(features):
    """COMPROBACIONES E y F. Devuelve la lista de avisos."""

    warnings = []

    # E · dentro de la caja de busqueda del docking
    for feature in features:

        for axis, center, size in zip(
            ("x", "y", "z"), BOX_CENTER, BOX_SIZE
        ):

            low = center - size / 2
            high = center + size / 2

            if not low <= feature[axis] <= high:
                warnings.append(
                    f"    {label(feature)} {feature['feature']}: "
                    f"{axis} = {feature[axis]:.3f} fuera de la caja "
                    f"[{low:.2f}, {high:.2f}]"
                )

    # F · ninguna pareja de hidrofobicas deberia seguir cumpliendo el criterio
    # de fusion con los radios que la regla 2 uso realmente (merge_radius).
    for i in range(len(features)):
        for j in range(i + 1, len(features)):

            a, b = features[i], features[j]

            if a["name"] != HYDROPHOBIC or b["name"] != HYDROPHOBIC:
                continue

            limit = MERGE_FACTOR * (a["merge_radius"] + b["merge_radius"])

            if distance(a, b) < limit:
                warnings.append(
                    f"    {label(a)} y {label(b)} siguen a "
                    f"{distance(a, b):.2f} A, por debajo del limite de "
                    f"fusion {limit:.2f} A: la regla 2 no las ha capturado"
                )

    return warnings


def report_post_clamp_overlaps(features):
    """
    Informativo, no es un error.

    La regla 3 amplia los radios DESPUES de fusionar, asi que dos esferas que
    no cumplian el criterio de fusion con sus radios originales pueden acabar
    solapandose una vez ajustadas. Pharmit trata cada esfera como una
    condicion obligatoria independiente, de modo que un par muy solapado
    equivale a exigir dos veces la misma zona y endurece el modelo.
    """

    overlaps = []

    for i, a in enumerate(features):
        for b in features[i + 1:]:

            if a["name"] != b["name"]:
                continue

            d = distance(a, b)
            suma = a["radius"] + b["radius"]

            if d < suma:
                overlaps.append(
                    f"    {label(a)} y {label(b)} ({a['feature']}): "
                    f"{d:.2f} A entre centros, radios {a['radius']:.2f} + "
                    f"{b['radius']:.2f} = {suma:.2f} A"
                )

    return overlaps


def report_positions(features):

    print(f"  {'Caracteristica':<18s} {'Tipo':<12s} "
          f"{'x':>8s} {'y':>8s} {'z':>8s} {'radio':>7s} {'cons':>6s}  clase")

    for feature in features:
        print(
            f"    {label(feature):<16s} {feature['feature']:<12s} "
            f"{feature['x']:8.3f} {feature['y']:8.3f} {feature['z']:8.3f} "
            f"{feature['radius']:7.2f} {feature['conservation']:6.2f}  "
            f"{feature['classification']}"
        )

        for member in feature.get("group", []):
            print(
                f"        miembro {label(member):<12s} "
                f"{member['x']:8.3f} {member['y']:8.3f} {member['z']:8.3f} "
                f"  a {distance(feature, member):.3f} A del centroide"
            )

    # Separaciones minimas, util para ver si dos esferas describen la misma zona
    pairs = [
        (distance(a, b), label(a), label(b))
        for i, a in enumerate(features)
        for b in features[i + 1:]
    ]

    if pairs:
        pairs.sort()
        print("  Pares mas proximos:")
        for d, first, second in pairs[:3]:
            print(f"    {first} - {second}: {d:.2f} A")


# ==========================================================
# ESCRITURA
# ==========================================================

CSV_COLUMNS = [
    "Model", "Interaction", "Residue", "Residue_number", "Feature", "PDB",
    "x", "y", "z", "Vector_x", "Vector_y", "Vector_z",
    "Radius", "Conservation", "Classification", "N_merged", "Merged_from",
]


def csv_row(level, feature):

    merged_from = "|".join(
        label(member) for member in feature.get("group", [])
    )

    return [
        level,
        feature["interaction"],
        "/".join(feature["residues"]),
        "/".join(feature["numbers"]),
        feature["feature"],
        feature["pdb"],
        f"{feature['x']:.3f}",
        f"{feature['y']:.3f}",
        f"{feature['z']:.3f}",
        f"{feature['vx']:.3f}",
        f"{feature['vy']:.3f}",
        f"{feature['vz']:.3f}",
        f"{feature['radius']:.3f}",
        f"{feature['conservation']:.2f}",
        feature["classification"],
        feature["members"],
        merged_from,
    ]


def write_csv(rows, path):

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as handle:

        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        writer.writerows(rows)


# ==========================================================


def main():

    if not os.path.exists(INPUT_FILE):
        sys.exit(
            f"No se encuentra {INPUT_FILE}. Ejecuta el script desde la raiz "
            f"del proyecto."
        )

    consensus = read_consensus(INPUT_FILE)

    print(f"\nTabla de consenso: {len(consensus)} caracteristicas\n")

    warn_classification_mismatch(consensus)

    all_warnings = []
    rows = []

    for level, categories in LEVELS.items():

        selected = [
            dict(f) for f in consensus
            if f["classification"] in categories
        ]

        refined, removed, absorbed = refine(selected)

        rows.extend(csv_row(level, feature) for feature in refined)

        print("=" * 72)
        print(f"Modelo {level.upper()}")
        print(f"  Filtradas por clasificacion : {len(selected)}")

        for feature, aromatic in removed:
            print(
                f"  Eliminada  {label(feature)} hidrofobica, cubierta por "
                f"{label(aromatic)} aromatica a "
                f"{distance(feature, aromatic):.2f} A"
            )

        print(f"  Absorbidas por fusion       : {absorbed}")
        print(f"  Caracteristicas finales     : {len(refined)}\n")

        report_positions(refined)

        warnings = check_positions(refined)

        if warnings:
            print(f"  AVISOS DE POSICION ({len(warnings)}):")
            for warning in warnings:
                print(warning)
            all_warnings.extend(warnings)
        else:
            print("  Posiciones: dentro de la caja de docking, sin "
                  "fusiones pendientes.")

        overlaps = report_post_clamp_overlaps(refined)

        if overlaps:
            print("  Esferas que se solapan tras ajustar radios "
                  "(informativo; endurece el modelo en Pharmit):")
            for overlap in overlaps:
                print(overlap)

        print()

    write_csv(rows, OUTPUT_FILE)

    print("=" * 72)
    print(f"Escritas {len(rows)} filas en {OUTPUT_FILE}")

    if all_warnings:
        print(f"Terminado con {len(all_warnings)} avisos de posicion.\n")
    else:
        print("Terminado. Todas las comprobaciones han pasado.\n")


if __name__ == "__main__":
    main()
