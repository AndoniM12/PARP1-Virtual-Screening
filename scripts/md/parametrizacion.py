from openff.toolkit import Molecule, ForceField
from openff.interchange import Interchange
from openff.units import unit
import os


# ============================================================
# CONFIGURACIÓN
# ============================================================

# Archivo SDF del ligando
ligand_file = (
    "data/md/ligand_453/ligando_453.sdf"
)

# Directorio de salida
# Salida a run/, que esta en .gitignore. La copia de la ejecucion
# original esta versionada en data/md/ligand_453/gromacs/.
output_dir = (
    "run/ligand/gromacs"
)

# Prefijo de los archivos GROMACS
output_prefix = os.path.join(
    output_dir,
    "ligand_453"
)

# Modelo NAGL para las cargas parciales
# Modelo NAGL para las cargas parciales.
#
# La ejecucion original apuntaba a una ruta absoluta dentro del entorno
# conda personal. Aqui se resuelve desde el paquete instalado para que
# no dependa de $HOME. Se puede forzar con la variable NAGL_MODEL.
#
NAGL_MODEL_NAME = "openff-gnn-am1bcc-1.0.0.pt"

nagl_model = os.environ.get("NAGL_MODEL")

if not nagl_model:
    try:
        from openff.nagl_models import get_nagl_model_path

        nagl_model = str(get_nagl_model_path(NAGL_MODEL_NAME))
    except Exception:
        nagl_model = NAGL_MODEL_NAME

print(f"Modelo NAGL: {nagl_model}")


# ============================================================
# CREAR DIRECTORIO DE SALIDA
# ============================================================

os.makedirs(output_dir, exist_ok=True)


# ============================================================
# CARGAR LIGANDO
# ============================================================

print("Cargando ligando...")

ligand = Molecule.from_file(
    ligand_file,
    file_format="SDF"
)

print(f"Átomos: {ligand.n_atoms}")
print(f"Nombre original: {ligand.name}")


# ============================================================
# ASIGNAR NOMBRE SIMPLE
# ============================================================
#
# El nombre original del SDF contiene múltiples identificadores
# separados por espacios. OpenFF utiliza este nombre durante la
# exportación a GROMACS, generando nombres de átomos/topología
# incompatibles con el formato esperado por GROMACS.
#
# Por ello utilizamos un nombre corto y estándar.
#

ligand.name = "LIG"

print(f"Nombre asignado: {ligand.name}")


# ============================================================
# CARGAR CAMPO DE FUERZA OPENFF
# ============================================================

print("\nCargando OpenFF 2.2.1...")

ff = ForceField(
    "openff-2.2.1.offxml"
)

print("Campo de fuerza cargado correctamente.")


# ============================================================
# CALCULAR CARGAS PARCIALES
# ============================================================

print("\nCalculando cargas parciales con NAGL...")

ligand.assign_partial_charges(
    partial_charge_method=nagl_model,
    normalize_partial_charges=True
)

print("Cargas parciales asignadas correctamente.")

print(
    "Carga total:",
    ligand.total_charge
)


# ============================================================
# CREAR TOPOLOGÍA OPENFF
# ============================================================

print("\nCreando topología del ligando...")

ligand_topology = ligand.to_topology()

print("Topología creada correctamente.")


# ============================================================
# CREAR INTERCHANGE
# ============================================================
#
# charge_from_molecules=[ligand] evita que Interchange vuelva
# a intentar calcular las cargas mediante otro método.
#

print("\nCreando sistema Interchange...")

interchange = Interchange.from_smirnoff(
    force_field=ff,
    topology=ligand_topology,
    charge_from_molecules=[ligand]
)

print("Interchange creado correctamente.")


# ============================================================
# CAJA TÉCNICA PARA LA EXPORTACIÓN
# ============================================================
#
# Estas dimensiones corresponden a la caja de docking:
#
# 22 x 16 x 24 Å
#
# convertidas a:
#
# 2.2 x 1.6 x 2.4 nm
#
# IMPORTANTE:
# Esta caja se utiliza únicamente para permitir la exportación
# a GROMACS. NO es la caja definitiva de la simulación MD.
#
# La caja definitiva se definirá posteriormente alrededor del
# complejo PARP1-ligando.
#

interchange.box = unit.Quantity(
    [
        [2.2, 0.0, 0.0],
        [0.0, 1.6, 0.0],
        [0.0, 0.0, 2.4],
    ],
    unit.nanometer
)

print("\nCaja técnica asignada:")
print("2.2 x 1.6 x 2.4 nm")


# ============================================================
# EXPORTAR A GROMACS
# ============================================================
#
# monolithic=False genera una topología modular, que es más
# adecuada para integrarla posteriormente con la topología
# de PARP1 generada mediante pdb2gmx.
#

print("\nExportando archivos GROMACS...")

interchange.to_gromacs(
    prefix=output_prefix,
    monolithic=False
)

print("\n============================================")
print("PARAMETRIZACIÓN COMPLETADA CORRECTAMENTE")
print("============================================")

print("\nArchivos generados en:")
print(output_dir)

print("\nPrefijo:")
print(output_prefix)