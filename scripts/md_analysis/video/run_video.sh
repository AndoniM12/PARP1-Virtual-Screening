#!/usr/bin/env bash
# Monta uno de los tres videos de principio a fin, sin interfaz grafica.
#
#   scripts/md_analysis/video/run_video.sh <animacion|amidas|torsiones>
#
# Desde la raiz del repositorio, y con data/md/traj_ligera.xtc en su sitio
# (se descarga del release; ver results/md/README.md).
#
# Existe porque la cadena son cuatro pasos repartidos entre ChimeraX y la
# terminal, y porque la forma de ejecutar ChimeraX sin interfaz --que es la
# que permite lanzarlo en segundo plano o en un servidor-- no estaba escrita
# en ninguna parte.
set -euo pipefail

VIDEO="${1:-}"
PYTHON="${PYTHON:-python}"
VD="scripts/md_analysis/video"
TRAJ="data/md/traj_ligera.xtc"

case "$VIDEO" in
  animacion)
    # La escena la monta el .cxc, que abre la estructura y suaviza por su cuenta.
    ESCENA="open $VD/estilo_revista.cxc"
    RENDER="$VD/render_animacion.py"
    PANEL="$VD/render_inset.py"
    COMPONER="$VD/componer.py"
    FINAL="run/render/final"
    FPS=60; CRF=16
    SALIDA="results/md/videos/ligando453_animacion_60fps.mp4"
    MINUTOS=20
    ;;
  amidas|torsiones)
    # Sin suaviza.py: la media movil amortiguaria las oscilaciones del diedro
    # que estos dos videos enseñan.
    ESCENA="open data/md/frame0.pdb; open $TRAJ structureModel #1"
    RENDER="$VD/render_$VIDEO.py"
    PANEL="$VD/inset_$VIDEO.py"
    COMPONER="$VD/componer_$VIDEO.py"
    FINAL="run/render/$VIDEO/final"
    SALIDA="results/md/videos/ligando453_$VIDEO.mp4"
    if [ "$VIDEO" = "amidas" ]; then FPS=30; CRF=18; MINUTOS=40;
    else FPS=60; CRF=18; MINUTOS=35; fi
    ;;
  *)
    echo "uso: $0 <animacion|amidas|torsiones>" >&2
    exit 2
    ;;
esac

for cmd in chimerax ffmpeg "$PYTHON"; do
  command -v "$cmd" >/dev/null || { echo "falta '$cmd' en el PATH" >&2; exit 1; }
done
[ -f "$TRAJ" ] || { echo "falta $TRAJ (se descarga del release)" >&2; exit 1; }
[ -d .git ] || { echo "ejecutalo desde la raiz del repositorio" >&2; exit 1; }

echo ">>> $VIDEO: renderizando con ChimeraX (~$MINUTOS min, deja ~500 MB en run/)"
LANZADOR="$(mktemp --suffix=.py)"
trap 'rm -f "$LANZADOR"' EXIT
{
  echo "from chimerax.core.commands import run"
  # El .cxc y los 'open' se pasan tal cual; el ; separa comandos de ChimeraX.
  echo "$ESCENA" | tr ';' '\n' | sed 's/^ *//' | while read -r c; do
    [ -n "$c" ] && printf 'run(session, %s)\n' "\"$c\""
  done
  printf 'run(session, "runscript %s")\n' "$RENDER"
} > "$LANZADOR"
# --silent se traga los session.logger.info de los scripts; los resumenes
# numericos salen de los scripts de analisis, no de estos.
chimerax --offscreen --nogui --silent --exit --script "$LANZADOR"

echo ">>> $VIDEO: panel"
"$PYTHON" "$PANEL"
echo ">>> $VIDEO: composicion"
"$PYTHON" "$COMPONER"

echo ">>> $VIDEO: codificando a $FPS fps"
mkdir -p "$(dirname "$SALIDA")"
ffmpeg -y -loglevel error -framerate "$FPS" -i "$FINAL/%04d.png" \
       -c:v libx264 -crf "$CRF" -pix_fmt yuv420p "$SALIDA"

echo
echo "listo: $SALIDA"
echo "run/ ocupa $(du -sh run 2>/dev/null | cut -f1); se puede borrar."
echo
echo "La miniatura del README se genera aparte, para no rehacerla sin querer:"
echo "  ffmpeg -i $SALIDA -vf 'fps=7,scale=720:-1:flags=lanczos' \\"
echo "         -c:v libwebp_anim -q:v 38 -loop 0 \\"
echo "         results/md/videos/ligando453_${VIDEO}_preview.webp"
