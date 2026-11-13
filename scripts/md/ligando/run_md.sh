#!/bin/bash
#
# Dinamica molecular del complejo PARP1 - un ligando cualquiera del cribado
# (50 ns, una replica). Es el protocolo del ligando 453 hecho repetible: el
# ligando se elige con LIGANDO, y lo unico que cambia entre dos ejecuciones
# es ese numero y la semilla del NVT.
#
# El 170 se simulo con este script; el 453 se hizo antes, a mano, y su
# ejecucion esta recuperada en scripts/md/run_md.sh. Los dos protocolos son
# el mismo (ver docs/02_molecular_dynamics.md), asi que una replica del 453
# lanzada desde aqui es comparable con la original.
#
# Uso, desde la raiz del repositorio:
#
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh prep
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh equil
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh prod
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh reanudar
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh estado
#     LIGANDO=170 bash scripts/md/ligando/run_md.sh todo     # equil y prod
#
# Para una tirada larga, lanzala de forma que sobreviva al cierre del
# terminal. Sin nohup, el proceso recibe SIGHUP y puede morir sin escribir
# el checkpoint del momento:
#
#     nohup env LIGANDO=170 bash scripts/md/ligando/run_md.sh todo \
#         > run/ligando_170/prod/run.out 2>&1 &
#     tail -f run/ligando_170/prod/run.out
#
# Hacen falta dos entornos de conda:
#   openff          para parametrizar el ligando (etapa prep)
#   md_parp1_cuda   para GROMACS (todo lo demas)
# El script los invoca con `conda run`, asi que no hay que activarlos a mano.
#
# Todo lo que se genera va a run/ligando_<N>/, que esta en .gitignore. Con
# -nb gpu la produccion de 50 ns tarda algo mas de un dia; en 4 hilos de CPU,
# unos cinco dias. El monitor termico para GROMACS limpiamente si la CPU pasa
# de TEMP_LIMIT grados; se reanuda con `reanudar` sin perder lo hecho.

set -euo pipefail

LIGANDO=${LIGANDO:-170}
STAGE=${1:-}

ENV_GMX=${ENV_GMX:-md_parp1_cuda}
ENV_OFF=${ENV_OFF:-openff}
CONDA=${CONDA:-"$HOME/miniconda3/bin/conda"}
NTOMP=${NTOMP:-4}
# Aceleracion de la produccion. El 170 corrio con "-nb gpu -pme cpu" (la
# GTX 1630 no compensa para PME) salvo los primeros 3,47 ns, que fueron en
# CPU. Con NB=cpu se reproduce la ejecucion del 453, que fue solo CPU.
NB=${NB:-gpu}
PME=${PME:-cpu}
TEMP_LIMIT=${TEMP_LIMIT:-80}
CHECK_INTERVAL=${CHECK_INTERVAL:-10}
# Cada cuantos MINUTOS se escribe el checkpoint. La simulacion del 453 uso 60,
# que es un corte de luz = una hora de calculo perdida. 15 es el valor por
# defecto de GROMACS y el coste de escribirlo es despreciable.
CPT=${CPT:-15}
# Semilla del generador de velocidades del NVT. Por defecto, el numero del
# ligando: asi queda escrita en el .mdp y la ejecucion es reproducible, al
# contrario que la del 453 (gen-seed = -1, ver docs/provenance.md). Con
# SEED=-1 se vuelve al comportamiento de aquella.
SEED=${SEED:-$LIGANDO}

WORK=${WORK:-run/ligando_$LIGANDO}
MDP="$WORK/mdp"
LIG="$WORK/ligand/ligando_$LIGANDO"
PROT="$WORK/protein"
CPLX="$WORK/complex"
EQ="$WORK"
PROD="$WORK/prod/md_50ns"

export WORK

gmx() { "$CONDA" run --no-capture-output -n "$ENV_GMX" gmx "$@"; }

mkdir -p "$WORK"/{mdp,protein,ligand,complex,minimization,nvt,npt,prod}

# Los .mdp son plantillas: @LIGANDO@ en el titulo y @SEED@ en la semilla del
# NVT. Se resuelven aqui, en $WORK/mdp, para que quede constancia de con que
# parametros exactos se lanzo cada ejecucion.
prepara_mdp() {
    local f
    for f in scripts/md/ligando/mdp/*.mdp; do
        sed -e "s/@LIGANDO@/$LIGANDO/g" -e "s/@SEED@/$SEED/g" \
            "$f" > "$MDP/$(basename "$f")"
    done
}

prep() {
    # La preparacion NO es repetible sobre si misma: solvate y genion escriben
    # en topol.top, asi que volver a lanzarla sobre una topologia ya completa
    # duplica lineas en [ molecules ]. Si hay que repetirla, se empieza limpio.
    if [ -f "$CPLX/complejo_iones.gro" ]; then
        echo "La preparacion ya esta hecha ($CPLX/complejo_iones.gro)." >&2
        echo "Para repetirla desde cero:" >&2
        echo "  rm -rf $PROT/* $CPLX/*   # el ligando parametrizado se conserva" >&2
        echo "Para seguir adelante:  LIGANDO=$LIGANDO bash scripts/md/ligando/run_md.sh equil" >&2
        exit 1
    fi
    prepara_mdp

    echo "== 1/5  Topologia de la proteina (amber19sb + TIP3P)"
    # 7AAA_repaired.pdb viene de la fase 1: solo cadena A, sin reconstruir los
    # residuos ausentes. Los hidrogenos los anade pdb2gmx.
    gmx pdb2gmx \
        -f data/md/topology/7AAA_repaired.pdb \
        -o "$PROT/PARP1.gro" \
        -p "$PROT/topol.top" \
        -i "$PROT/posre.itp" \
        -ff amber19sb \
        -water tip3p

    echo
    echo "== 2/5  Ligando en la pose de GNINA (OpenFF 2.2.1 + cargas NAGL)"
    if [ -f "$LIG.gro" ]; then
        echo "   ya esta hecho ($LIG.gro); se reutiliza"
    else
        "$CONDA" run --no-capture-output -n "$ENV_OFF" \
            python scripts/md/ligando/preparar_ligando.py "$LIGANDO"
    fi

    echo
    echo "== 3/5  Montar el complejo y meter el ligando en la topologia"
    "$CONDA" run --no-capture-output -n "$ENV_OFF" \
        python scripts/md/ligando/montar_complejo.py "$LIGANDO"

    echo
    echo "== 4/5  Caja dodecaedrica (1,0 nm de margen) y solvatacion"
    gmx editconf \
        -f "$CPLX/complejo.gro" \
        -o "$CPLX/complejo_caja.gro" \
        -bt dodecahedron \
        -d 1.0

    local spc
    spc="$("$CONDA" run -n "$ENV_GMX" \
           bash -c 'echo "$CONDA_PREFIX/share/gromacs/top/spc216.gro"')"
    gmx solvate \
        -cp "$CPLX/complejo_caja.gro" \
        -cs "$spc" \
        -o "$CPLX/complejo_solvatado.gro" \
        -p "$PROT/topol.top"

    echo
    echo "== 5/5  Neutralizar con contraiones"
    # maxwarn 1: el aviso esperado por la carga neta antes de anadir iones.
    gmx grompp \
        -f "$MDP/ions.mdp" \
        -c "$CPLX/complejo_solvatado.gro" \
        -p "$PROT/topol.top" \
        -o "$CPLX/ions.tpr" \
        -maxwarn 1
    # SOL: el grupo del que se sustituyen moleculas por iones.
    echo SOL | gmx genion \
        -s "$CPLX/ions.tpr" \
        -o "$CPLX/complejo_iones.gro" \
        -p "$PROT/topol.top" \
        -pname NA -nname CL -neutral

    echo
    echo "PREPARACION TERMINADA. Siguiente:"
    echo "  LIGANDO=$LIGANDO bash scripts/md/ligando/run_md.sh equil"
}

equil() {
    prepara_mdp
    echo "== Minimizacion de energia (steep, emtol 1000)"
    gmx grompp -f "$MDP/em.mdp" -c "$CPLX/complejo_iones.gro" \
        -p "$PROT/topol.top" -o "$EQ/minimization/em.tpr"
    gmx mdrun -deffnm "$EQ/minimization/em" -ntmpi 1 -ntomp "$NTOMP" -nb cpu -v

    echo
    echo "== Equilibrado NVT (100 ps, 300 K, semilla $SEED)"
    gmx grompp -f "$MDP/nvt.mdp" -c "$EQ/minimization/em.gro" \
        -p "$PROT/topol.top" -o "$EQ/nvt/nvt.tpr"
    gmx mdrun -deffnm "$EQ/nvt/nvt" -ntmpi 1 -ntomp "$NTOMP" -nb cpu -v

    echo
    echo "== Equilibrado NPT (100 ps, 1 bar)"
    gmx grompp -f "$MDP/npt.mdp" -c "$EQ/nvt/nvt.gro" \
        -p "$PROT/topol.top" -o "$EQ/npt/npt.tpr"
    gmx mdrun -deffnm "$EQ/npt/npt" -ntmpi 1 -ntomp "$NTOMP" -nb cpu -v

    echo
    echo "EQUILIBRADO TERMINADO."
    echo "Los .xtc de NVT y NPT se CONSERVAN a proposito: de la simulacion del"
    echo "453 no se guardaron y por eso no se pudo fechar en que momento se"
    echo "perdio el puente de hidrogeno con Gly863. No los borres."
    echo
    echo "Siguiente:  LIGANDO=$LIGANDO bash scripts/md/ligando/run_md.sh prod"
}

prod() {
    prepara_mdp
    if [ ! -f "$PROD.tpr" ]; then
        echo "== Preparando la produccion de 50 ns"
        gmx grompp \
            -f "$MDP/md_50ns.mdp" \
            -c "$EQ/npt/npt.gro" \
            -r "$EQ/npt/npt.gro" \
            -t "$EQ/npt/npt.cpt" \
            -p "$PROT/topol.top" \
            -o "$PROD.tpr"
    fi
    lanzar_con_monitor
}

reanudar() {
    [ -f "$PROD.cpt" ] || { echo "no hay $PROD.cpt que reanudar" >&2; exit 1; }
    echo "== Reanudando desde el checkpoint"
    EXTRA=(-cpi "$PROD.cpt" -append)
    lanzar_con_monitor
}

lanzar_con_monitor() {
    echo "Monitor termico activo. Limite ${TEMP_LIMIT} C, comprobando cada ${CHECK_INTERVAL} s."
    gmx mdrun \
        -s "$PROD.tpr" \
        -deffnm "$PROD" \
        ${EXTRA[@]+"${EXTRA[@]}"} \
        -ntmpi 1 -ntomp "$NTOMP" -nb "$NB" -pme "$PME" -cpt "$CPT" -v &
    local pid=$!
    echo "GROMACS lanzado con PID $pid"

    while kill -0 "$pid" 2>/dev/null; do
        local t
        t=$(sensors 2>/dev/null | grep -E 'Core [0-9]+:' \
            | awk '{gsub(/[+°C]/,"",$3); print $3}' | sort -nr | head -1)
        if [ -n "${t:-}" ]; then
            printf '%s - CPU maxima: %s C\n' "$(date '+%H:%M:%S')" "$t"
            if awk "BEGIN {exit !($t >= $TEMP_LIMIT)}"; then
                echo "TEMPERATURA ALTA: ${t} C >= ${TEMP_LIMIT} C. Parando limpiamente."
                kill -INT "$pid"
                wait "$pid" || true
                echo "Detenido. Continuar con:"
                echo "  LIGANDO=$LIGANDO bash scripts/md/ligando/run_md.sh reanudar"
                return 0
            fi
        fi
        sleep "$CHECK_INTERVAL"
    done
    wait "$pid" || true
    echo "PRODUCCION TERMINADA."
}

estado() {
    echo "ligando $LIGANDO, directorio de trabajo $WORK"
    printf '%-34s %s\n' "ligando parametrizado" \
        "$([ -f "$LIG.gro" ] && echo si || echo NO)"
    printf '%-34s %s\n' "topologia de la proteina" \
        "$([ -f "$PROT/topol.top" ] && echo si || echo NO)"
    printf '%-34s %s\n' "complejo solvatado y neutro" \
        "$([ -f "$CPLX/complejo_iones.gro" ] && echo si || echo NO)"
    local e d n
    for e in minimization:em nvt:nvt npt:npt; do
        d=${e%%:*}; n=${e##*:}
        printf '%-34s %s\n' "equilibrado: $d" \
            "$([ -f "$EQ/$d/$n.gro" ] && echo si || echo NO)"
    done
    printf '%-34s %s\n' "produccion lanzada" \
        "$([ -f "$PROD.log" ] && echo si || echo NO)"
    if [ -f "$PROD.log" ]; then
        echo
        echo "ultimo paso registrado:"
        grep -E "^ *Step +Time" -A1 "$PROD.log" | tail -2
        echo "(50 ns = 25 000 000 pasos)"
    fi
}

case "$STAGE" in
    prep)     prep ;;
    equil)    equil ;;
    prod)     prod ;;
    reanudar) reanudar ;;
    todo)     equil; prod ;;
    estado)   estado ;;
    *)
        sed -n '13,24p' "$0" | sed 's/^# \{0,1\}//'
        exit 1 ;;
esac
