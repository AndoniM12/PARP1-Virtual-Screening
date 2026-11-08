# Entorno del analisis de la dinamica (ChimeraX)

ChimeraX no se instala con conda, asi que no tiene un `.yml`.

## Version

UCSF ChimeraX **1.12**, instalado desde el `.deb` oficial para
Ubuntu 24.04 (`ucsf-chimerax_1.12ubuntu24.04_amd64.deb`).

La version se ha deducido del nombre del paquete instalado que consta
en el historial de terminal; no estaba declarada en ningun fichero del
proyecto.

## Scripts que lo necesitan

Se ejecutan **dentro de ChimeraX**, con `runscript`:

- `scripts/md_analysis/analiza.py` (usa `chimerax.hbonds.find_hbonds`)
- `scripts/md_analysis/anillos.py`
- `scripts/md_analysis/suaviza.py`
- `scripts/md_analysis/render_animacion.py`
- `scripts/md_analysis/estilo_revista.cxc`

## Scripts que NO lo necesitan

Son Python normal y corriente, con matplotlib, numpy y Pillow. Valen
tanto el entorno `parp1_pipeline` como cualquiera que tenga esos tres
paquetes:

- `scripts/md_analysis/grafica.py`
- `scripts/md_analysis/render_inset.py`
- `scripts/md_analysis/componer.py`

## Otras dependencias

- `ffmpeg`, para montar el video a partir de los fotogramas.
- Fuentes DejaVu en `/usr/share/fonts/truetype/dejavu/`, que
  `componer.py` usa por ruta absoluta.

## Aviso sobre las rutas

ChimeraX resuelve las rutas relativas respecto a su directorio de
trabajo. Arrancalo desde la raiz del repositorio, o ejecuta antes el
comando `cd <raiz del repositorio>` en su linea de comandos.
