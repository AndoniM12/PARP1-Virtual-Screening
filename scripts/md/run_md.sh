#!/bin/bash
#
# Dinamica molecular del complejo PARP1 - ligando 453 (50 ns, una replica).
#
# Reproduce la ejecucion original. Los pasos 1 a 9 se han recuperado del
# historial de terminal: no existian como script. Las rutas son relativas
# a la raiz del repositorio.
#
# Uso:
#     bash scripts/md/run_md.sh prep      # pasos 1-9 (preparacion)
#     bash scripts/md/run_md.sh equil     # pasos 10-12 (minimizacion y equilibrado)
#     bash scripts/md/run_md.sh prod      # paso 13 (produccion de 50 ns)
#     bash scripts/md/run_md.sh analysis  # paso 14 (centrado de la trayectoria)
#
# Todo lo que se genera va a run/, que esta en .gitignore. La produccion
# tardo unos 5 dias en 4 hilos de CPU (ver docs/02_molecular_dynamics.md).
#
# AVISO: la produccion original se lanzo con semilla aleatoria
# (gen-seed = -1). Para reproducirla exactamente hay que fijar las
# semillas efectivas, que estan en docs/provenance.md.

set -euo pipefail

WORK=${WORK:-run}
STAGE=${1:-}

# Limite termico del monitor de la produccion (grados C).
TEMP_LIMIT=${TEMP_LIMIT:-80}
CHECK_INTERVAL=${CHECK_INTERVAL:-10}

# Hilos. La ejecucion original uso 4 hilos OpenMP y solo CPU.
NTOMP=${NTOMP:-4}

mkdir -p "$WORK"/{protein,ligand,complex,minimization,nvt,npt,prod}

prep() {
    echo "== 1. Reparar la estructura experimental (PDBFixer)"
    # Solo cadena A. No se reconstruyen residuos ausentes (Gly660/Ser661
    # quedan fuera a proposito). Los H los anade pdb2gmx despues.
    python scripts/md/reparar_proteina.py

    echo "== 2. Topologia de la proteina (amber19sb + TIP3P)"
    gmx pdb2gmx \
        -f data/md/topology/7AAA_repaired.pdb \
        -o "$WORK/protein/PARP1.gro" \
        -p "$WORK/protein/topol.top" \
        -i "$WORK/protein/posre.itp" \
        -ff amber19sb \
        -water tip3p

    echo "== 3. Parametrizar el ligando (OpenFF 2.2.1 + cargas NAGL)"
    # Requiere el entorno `openff`, no `md_parp1`.
    # Ver environment/ligand_openff.yml
    python scripts/md/parametrizacion.py

    echo "== 4. Comprobar que el SDF conserva la pose del docking"
    python scripts/md/comprobar_pose_ligando.py

    echo "== 5. Mapear SDF <-> MOL2 de GNINA"
    python scripts/md/construir_ligando_docked.py

    echo "== 6. Generar el .gro del ligando en la pose dockeada"
    python scripts/md/generar_ligando_docked_gro.py

    echo "== 7. Validar la geometria de los hidrogenos"
    python scripts/md/validar_ligando_docked.py

    echo "== 8. Anadir el ligando a [ molecules ] de la topologia"
    # En la ejecucion original esto se hizo con sed sobre topol.top.
    sed -i '/Protein_chain_A     1/a LIG                 1' "$WORK/protein/topol.top"

    echo "== 9. Montar el complejo, caja, solvatacion e iones"
    python scripts/md/montar_complejo.py "$WORK"

    gmx editconf \
        -f "$WORK/complex/complex.gro" \
        -o "$WORK/complex/complex_box.gro" \
        -bt dodecahedron \
        -d 1.0

    gmx solvate \
        -cp "$WORK/complex/complex_box.gro" \
        -cs "$CONDA_PREFIX/share/gromacs/top/spc216.gro" \
        -o "$WORK/complex/complex_solvated.gro" \
        -p "$WORK/protein/topol.top"

    # maxwarn 1: aviso esperado por la carga neta del sistema antes de
    # anadir contraiones.
    gmx grompp \
        -f scripts/md/mdp/ions.mdp \
        -c "$WORK/complex/complex_solvated.gro" \
        -p "$WORK/protein/topol.top" \
        -o "$WORK/complex/ions.tpr" \
        -maxwarn 1

    # Neutraliza el sistema. En la ejecucion original resulto 1 solo Cl-.
    # No se anadio sal a fuerza ionica fisiologica (ver limitaciones).
    gmx genion \
        -s "$WORK/complex/ions.tpr" \
        -o "$WORK/complex/complex_ions.gro" \
        -p "$WORK/protein/topol.top" \
        -pname NA \
        -nname CL \
        -neutral
}

equil() {
    echo "== 10. Minimizacion de energia"
    gmx grompp \
        -f scripts/md/mdp/em.mdp \
        -c "$WORK/complex/complex_ions.gro" \
        -p "$WORK/protein/topol.top" \
        -o "$WORK/minimization/em.tpr"
    gmx mdrun -deffnm "$WORK/minimization/em" -nt "$NTOMP"

    echo "== 11. Equilibrado NVT (100 ps, 300 K)"
    gmx grompp \
        -f scripts/md/mdp/nvt.mdp \
        -c "$WORK/minimization/em.gro" \
        -p "$WORK/protein/topol.top" \
        -o "$WORK/nvt/nvt.tpr"
    gmx mdrun -deffnm "$WORK/nvt/nvt" -nt "$NTOMP"

    echo "== 12. Equilibrado NPT (100 ps, 1 bar)"
    gmx grompp \
        -f scripts/md/mdp/npt.mdp \
        -c "$WORK/nvt/nvt.gro" \
        -p "$WORK/protein/topol.top" \
        -o "$WORK/npt/npt.tpr"
    gmx mdrun -deffnm "$WORK/npt/npt" -nt "$NTOMP"
}

prod() {
    echo "== 13. Produccion, 50 ns"
    gmx grompp \
        -f scripts/md/mdp/md_50ns.mdp \
        -c "$WORK/npt/npt.gro" \
        -r "$WORK/npt/npt.gro" \
        -t "$WORK/npt/npt.cpt" \
        -p "$WORK/protein/topol.top" \
        -o "$WORK/prod/md_50ns.tpr"

    local prefix="$WORK/prod/md_50ns"

    # Monitor termico: la ejecucion original se hizo en un portatil sin
    # GPU y se paraba limpiamente si la CPU superaba el limite. GROMACS
    # se reanuda despues con `gmx mdrun -cpi`.
    echo "Monitor termico activo. Limite: ${TEMP_LIMIT} C"

    gmx mdrun \
        -s "${prefix}.tpr" \
        -deffnm "${prefix}" \
        -ntmpi 1 \
        -ntomp "$NTOMP" \
        -nb cpu \
        -cpt 60 \
        -v &

    local gmx_pid=$!
    echo "GROMACS lanzado con PID $gmx_pid"

    while kill -0 "$gmx_pid" 2>/dev/null; do
        local temp
        temp=$(sensors 2>/dev/null | grep -E 'Core [0-9]+:' \
            | awk '{gsub(/[+°C]/,"",$3); print $3}' | sort -nr | head -1)

        if [ -n "${temp:-}" ]; then
            echo "$(date '+%H:%M:%S') - CPU maxima: ${temp} C"
            if awk "BEGIN {exit !($temp >= $TEMP_LIMIT)}"; then
                echo "TEMPERATURA ALTA: ${temp} C >= ${TEMP_LIMIT} C. Parando."
                kill -INT "$gmx_pid"
                wait "$gmx_pid" || true
                echo "Detenido. Reanudar con: gmx mdrun -s ${prefix}.tpr -cpi ${prefix}.cpt -deffnm ${prefix}"
                return 0
            fi
        fi
        sleep "$CHECK_INTERVAL"
    done

    wait "$gmx_pid" || true
    echo "Produccion terminada."
}

analysis() {
    echo "== 14. Centrar la trayectoria"
    # Grupo a elegir de forma interactiva en la ejecucion original.
    gmx trjconv \
        -s "$WORK/prod/md_50ns.tpr" \
        -f "$WORK/prod/md_50ns.xtc" \
        -o "$WORK/prod/md_50ns_centered.xtc" \
        -pbc mol \
        -center

    echo
    echo "NO RECUPERADO: el comando que genero data/md/ligand_453/traj_ligera.xtc"
    echo "(sin agua, ajustada al backbone, 1 fotograma cada 50 ps) no"
    echo "aparece en el historial. Ver docs/provenance.md."
}

case "$STAGE" in
    prep)     prep ;;
    equil)    equil ;;
    prod)     prod ;;
    analysis) analysis ;;
    *)
        echo "Uso: bash scripts/md/run_md.sh {prep|equil|prod|analysis}"
        exit 1
        ;;
esac
