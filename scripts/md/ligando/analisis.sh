#!/bin/bash
#
# Analisis con GROMACS de una produccion de 50 ns: centra la trayectoria,
# la alinea, y saca RMSD, RMSF, radio de giro, distancia al bolsillo,
# puentes de hidrogeno y las energias de control.
#
# Es lo que se ejecuto para el ligando 170. Las tablas van a
# results/md/ligand_<N>/tables/ y la trayectoria ligera a
# data/md/ligand_<N>/traj_fit.xtc (que esta en .gitignore y se publica
# como release).
#
# Uso, desde la raiz del repositorio:
#
#     LIGANDO=170 bash scripts/md/ligando/analisis.sh
#
# Necesita la produccion en run/ligando_<N>/prod/ (el .tpr y el .xtc
# completos, que no estan en el repositorio por tamano). Con SRC se puede
# apuntar a otro sitio:
#
#     SRC=/ruta/a/05_produccion LIGANDO=170 bash scripts/md/ligando/analisis.sh
#
# El indice se construye aqui mismo: ademas de los grupos por defecto hace
# falta Proteina_LIG (proteina + ligando, para centrar y escribir sin agua)
# y LIG_pesados (los atomos pesados del ligando, que son los que tienen
# coordenadas del docking y con los que se mide el RMSD).

set -euo pipefail

LIGANDO=${LIGANDO:-170}
WORK=${WORK:-run/ligando_$LIGANDO}
SRC=${SRC:-$WORK/prod}
PROD=${PROD:-$SRC/md_50ns}
OUT=${OUT:-results/md/ligand_$LIGANDO/tables}
DATA=${DATA:-data/md/ligand_$LIGANDO}
# Paso entre fotogramas de la trayectoria ligera, en ps. 100 deja 501
# fotogramas de 50 ns, que es lo que se uso con el 170.
DT=${DT:-100}
# Residuos del bolsillo de la nicotinamida con los que se mide la distancia.
BOLSILLO=${BOLSILLO:-"863 904 907"}

ENV_GMX=${ENV_GMX:-md_parp1_cuda}
CONDA_PREFIX_GMX=${CONDA_PREFIX_GMX:-"$HOME/miniconda3/envs/$ENV_GMX"}
# Se invoca el binario del entorno directamente, no con `conda run`: este
# no reenvia la entrada estandar, y casi todas las herramientas de aqui
# piden los grupos por ahi.
GMX="$CONDA_PREFIX_GMX/bin/gmx -quiet"

for f in "$PROD.tpr" "$PROD.xtc" "$PROD.edr"; do
    [ -f "$f" ] || { echo "falta $f (define SRC o PROD)" >&2; exit 1; }
done
mkdir -p "$OUT/equilibration" "$DATA" "$WORK/logs"
L="$WORK/logs"

echo "== Indice: Proteina_LIG y LIG_pesados"
printf '1 | 13\nname 20 Proteina_LIG\n13 & ! a H*\nname 21 LIG_pesados\nq\n' \
    | $GMX make_ndx -f "$PROD.tpr" -o "$WORK/index.ndx" > "$L/mk.log" 2>&1
NDX="$WORK/index.ndx"
# Numeros de grupo, ya con los dos anadidos.
G_PROT=1; G_BB=4; G_CA=3; G_LIG=13; G_PLIG=20; G_LIGP=21

echo "== Centrar en la proteina y escribir sin agua (1 fotograma cada $DT ps)"
printf '%s\n%s\n' $G_PROT $G_PLIG | $GMX trjconv -s "$PROD.tpr" -f "$PROD.xtc" \
    -n "$NDX" -o "$WORK/centrada.xtc" -pbc mol -center -ur compact -dt "$DT" \
    > "$L/trj.log" 2>&1

echo "== Alinear por el backbone"
printf '%s\n%s\n' $G_BB $G_PLIG | $GMX trjconv -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$DATA/traj_fit.xtc" -fit rot+trans > "$L/fit.log" 2>&1
printf '%s\n%s\n' $G_BB $G_PLIG | $GMX trjconv -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$DATA/frame0.pdb" -fit rot+trans -dump 0 > "$L/pdb.log" 2>&1

# El RMSD se mide sobre la trayectoria centrada, no sobre la alineada: gmx
# rms hace su propio ajuste con el grupo de minimos cuadrados.
echo "== RMSD"
printf '%s\n%s\n' $G_BB $G_BB | $GMX rms -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsd_backbone.xvg" -tu ns > "$L/rms_bb.log" 2>&1
# Ligando con la proteina alineada: cuanto se ha movido en el bolsillo.
printf '%s\n%s\n' $G_BB $G_LIGP | $GMX rms -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsd_ligpes_en_bolsillo.xvg" -tu ns > "$L/rms_ligp.log" 2>&1
# Ligando sobre si mismo: solo la deformacion interna.
printf '%s\n%s\n' $G_LIGP $G_LIGP | $GMX rms -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsd_ligpes_propio.xvg" -tu ns > "$L/rms_ligp_self.log" 2>&1
# Los mismos dos, con los hidrogenos incluidos.
printf '%s\n%s\n' $G_BB $G_LIG | $GMX rms -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsd_lig_en_bolsillo.xvg" -tu ns > "$L/rms_lig.log" 2>&1
printf '%s\n%s\n' $G_LIG $G_LIG | $GMX rms -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsd_lig_propio.xvg" -tu ns > "$L/rms_lig_self.log" 2>&1

echo "== RMSF por residuo (5-50 ns, saltando el acomodo inicial)"
printf '%s\n' $G_CA | $GMX rmsf -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rmsf_ca.xvg" -res -b 5000 > "$L/rmsf.log" 2>&1

echo "== Radio de giro"
printf '%s\n' $G_PROT | $GMX gyrate -s "$PROD.tpr" -f "$WORK/centrada.xtc" \
    -n "$NDX" -o "$OUT/rg.xvg" > "$L/rg.log" 2>&1

echo "== Distancia del ligando al bolsillo"
$GMX distance -s "$PROD.tpr" -f "$WORK/centrada.xtc" -n "$NDX" \
    -select "com of resname LIG plus com of (resid $BOLSILLO and name CA)" \
    -oall "$OUT/dist_bolsillo.xvg" -tu ns > "$L/dist.log" 2>&1

echo "== Puentes de hidrogeno ligando-proteina"
# Criterio de gmx hbond (0,35 nm y 30 grados). No coincide con el de
# ChimeraX, que es el que usan las tablas de scripts/md_analysis.
$GMX hbond -s "$PROD.tpr" -f "$WORK/centrada.xtc" -n "$NDX" \
    -r LIG -t Protein -num "$OUT/hbnum.xvg" -tu ns > "$L/hb.log" 2>&1 \
    || echo "   gmx hbond ha fallado; sigue el resto" >&2

echo "== Energias de control de la produccion"
for t in Temperature Pressure Density Potential; do
    echo "$t" | $GMX energy -f "$PROD.edr" -o "$OUT/ener_$t.xvg" \
        > "$L/ener_$t.log" 2>&1
done

echo "== Equilibrado"
if [ -f "$WORK/nvt/nvt.edr" ]; then
    echo Temperature | $GMX energy -f "$WORK/nvt/nvt.edr" \
        -o "$OUT/equilibration/nvt_temperature.xvg" > "$L/eq_nvt.log" 2>&1
fi
if [ -f "$WORK/npt/npt.edr" ]; then
    echo Temperature | $GMX energy -f "$WORK/npt/npt.edr" \
        -o "$OUT/equilibration/npt_temperature.xvg" > "$L/eq_npt_t.log" 2>&1
    printf 'Pressure\nDensity\n' | $GMX energy -f "$WORK/npt/npt.edr" \
        -o "$OUT/equilibration/npt_pressure_density.xvg" > "$L/eq_npt_pd.log" 2>&1
    printf '%s\n%s\n' $G_BB $G_BB | $GMX rms -s "$WORK/npt/npt.tpr" \
        -f "$WORK/npt/npt.xtc" -n "$NDX" \
        -o "$OUT/equilibration/npt_rmsd_backbone.xvg" > "$L/eq_npt_rms.log" 2>&1
fi

echo
echo "TERMINADO. Tablas en $OUT, trayectoria ligera en $DATA/traj_fit.xtc"
echo "Siguiente, en ChimeraX desde la raiz del repositorio:"
echo "  open $DATA/frame0.pdb"
echo "  open $DATA/traj_fit.xtc structureModel #1"
echo "  runscript scripts/md_analysis/analiza.py $LIGANDO 0.1"
