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
- `scripts/md_analysis/amidas.py`
- `scripts/md_analysis/torsiones.py`
- `scripts/md_analysis/fragmentos.py`
- `scripts/md_analysis/conformeros.py`

## Scripts que NO lo necesitan

Son Python normal y corriente, con numpy:

- `scripts/md_analysis/escena_cuatro.py`
- `scripts/md_analysis/tension.py`
- `scripts/md_analysis/alcance.py`

`escena_cuatro.py` escribe una escena `.cxc` que sí se abre con
ChimeraX, pero el script en sí no lo necesita para correr.

## Aviso sobre las rutas

ChimeraX resuelve las rutas relativas respecto a su directorio de
trabajo. Arrancalo desde la raiz del repositorio, o ejecuta antes el
comando `cd <raiz del repositorio>` en su linea de comandos.

## Que ligando analizan

Los cinco scripts que leen una trayectoria (`analiza.py`, `anillos.py`,
`amidas.py`, `torsiones.py` y `fragmentos.py`) toman dos argumentos: el
numero del ligando y el paso entre fotogramas en ns. Sin argumentos
analizan el 453 con sus 50 ps, que es como estaban antes.

    runscript scripts/md_analysis/analiza.py          # 453, 50 ps
    runscript scripts/md_analysis/analiza.py 170 0.1  # 170, 100 ps

Las tablas salen a `results/md/ligand_<N>/tables/`. El ligando NO se lee
de la estructura cargada: los scripts analizan lo que haya en `#1`, asi
que si los argumentos no se corresponden con la trayectoria abierta, las
tablas salen con el nombre equivocado y la escala de tiempo mal.

`amidas.py`, `torsiones.py` y `fragmentos.py` identifican los grupos por
la quimica del **453** (sus dos amidas, su enlazador, su benzotriazol),
asi que con otro ligando fallan o no encuentran nada. Con el 170 solo
valen `analiza.py` (puentes de H y contactos) y `anillos.py`.
