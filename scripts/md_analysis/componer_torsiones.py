# Compone cada fotograma del video de las torsiones: render de ChimeraX a la
# izquierda, panel de diales a la derecha, titulo y reloj.
# Salida: run/render/torsiones/final/NNNN.png
#
#   python scripts/md_analysis/componer_torsiones.py
import csv
import os
import sys
import textwrap

from PIL import Image, ImageDraw, ImageFont

# Por defecto, la raiz del repositorio.
ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
D = ROOT + "/run/render/torsiones/"        # en .gitignore
os.makedirs(D + "final", exist_ok=True)

ANCHO, ALTO = 2000, 900                    # dimensiones pares: ffmpeg las exige
FUENTE = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
NEGRITA = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
titulo = ImageFont.truetype(NEGRITA, 30)
reloj = ImageFont.truetype(FUENTE, 26)
pie = ImageFont.truetype(FUENTE, 18)
aviso = ImageFont.truetype(NEGRITA, 21)

NOTA = textwrap.fill(
    "El núcleo (gris) está fijo: no cambia de forma en 50 ns. Lo que se mueve "
    "es el brazo del benzotriazol (violeta claro), y lo mueven los tres "
    "enlaces coloreados. El decisivo es el CH₂–CH₂: salta en 27,75 y 30,65 ns "
    "y manda el brazo de 2 a 7 Å.", width=62)

for p in csv.DictReader(open(D + "plan.csv")):
    n = int(p["salida"])
    t = float(p["tiempo_ns"])
    mol = Image.open(f"{D}mol/{n:04d}.png").convert("RGBA")
    inset = Image.open(f"{D}inset/{n:04d}.png").convert("RGBA")

    img = Image.new("RGBA", (ANCHO, ALTO), (255, 255, 255, 255))
    img.alpha_composite(mol, (0, (ALTO - mol.height) // 2))
    x0 = mol.width + (ANCHO - mol.width - inset.width) // 2
    img.alpha_composite(inset, (x0, 150))

    d = ImageDraw.Draw(img)
    d.line((mol.width + 10, 70, mol.width + 10, ALTO - 70), fill="#e3e8ee", width=2)
    d.text((44, 36), "PARP1 · ligando 453 · qué mueve el cambio de modo de unión",
           font=titulo, fill="#1f2933")
    d.text((44, 80), f"t = {t:.1f} ns", font=reloj, fill="#52606d")
    d.text((x0, ALTO - 185), NOTA, font=pie, fill="#52606d", spacing=7)
    if 27.5 <= t <= 30.0:
        d.rounded_rectangle((40, 124, 560, 162), radius=8, fill="#ede9fe")
        d.text((52, 130), "cambio de modo de unión", font=aviso, fill="#6d28d9")
    img.convert("RGB").save(f"{D}final/{n:04d}.png")
print(f"Compuestos en {D}final/")
