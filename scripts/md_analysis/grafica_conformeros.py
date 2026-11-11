# Compone la figura de las cuatro etapas conformacionales del ligando 453 a
# partir de las parejas que renderiza conformeros.py.
# Salida: results/md/figures/conformeros.png
#
#   python scripts/md_analysis/grafica_conformeros.py
import os
import sys
import textwrap

from PIL import Image, ImageDraw, ImageFont

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
RENDER = ROOT + "/run/render/conformeros/"      # en .gitignore
FIGURAS = ROOT + "/results/md/figures"
os.makedirs(FIGURAS, exist_ok=True)

FUENTE = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
NEGRITA = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
titulo = ImageFont.truetype(NEGRITA, 30)
leyenda = ImageFont.truetype(FUENTE, 22)
pie = ImageFont.truetype(FUENTE, 21)

PANELES = [
    ("par_1_docking.png", "1 · El docking retuerce el enlazador",
     [("biblioteca de Pharmit", "#2563eb", "CH₂–CH₂ = −179,9°  (anti)"),
      ("pose de GNINA", "#ea580c", "CH₂–CH₂ = 115,2°  (casi eclipsado)")]),
    ("par_2_minimizacion.png", "2 · La minimización alivia el eclipse",
     [("pose de GNINA", "#ea580c", "115,2°"),
      ("tras minimizar y equilibrar", "#0f766e", "77,4°  ·  cae al mínimo vecino")]),
    ("par_3_dinamica.png", "3 · La dinámica cruza a anti",
     [("inicio de la dinámica", "#0f766e", "77,4°  (g+)"),
      ("t = 50 ns", "#7c3aed", "−171,0°  (anti)")]),
    ("par_4_vuelta.png", "4 · El diedro vuelve, el brazo no",
     [("biblioteca de Pharmit", "#2563eb", "−179,9°"),
      ("t = 50 ns", "#7c3aed", "−171,0°  ·  pero CO–CH₂ acabó en −8,5°")])]

# Los valores son los de las estructuras que se dibujan, instantaneos. Una
# version anterior etiquetaba aqui las medias circulares de las dos mesetas
# (64,2° y -172,1°), que no son el valor de ningun fotograma concreto: el panel
# 2 y el 3 ensenan la misma estructura y la daban con dos numeros distintos.

NOTA = textwrap.fill(
    "Las cuatro conformaciones superpuestas sobre el anillo central: RMSD del "
    "núcleo ≤ 0,10 Å, y 0,000 Å en la pose de GNINA, que mantiene rígido el "
    "esqueleto y solo mueve torsiones. Ligando 453 · MolPort-000-056-093.",
    width=132)

ancho, alto = 1400, 1000
esc = 0.62
cw, ch = int(ancho * esc), int(alto * esc)
alto_leyenda = 2 * 34 + 14
W = cw * 2 + 60
H = (ch + alto_leyenda + 60) * 2 + 80

img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)


def ancho_texto(txt, fuente):
    """La anchura real del texto. Fijar una columna a ojo hacia que las
    etiquetas largas se comieran los valores."""
    return d.textbbox((0, 0), txt, font=fuente)[2]


for k, (fichero, tit, leyendas) in enumerate(PANELES):
    x = 20 + (k % 2) * (cw + 20)
    y = 60 + (k // 2) * (ch + alto_leyenda + 60)
    img.paste(Image.open(RENDER + fichero).convert("RGB")
              .resize((cw, ch), Image.LANCZOS), (x, y))
    d.text((x + 6, y - 42), tit, font=titulo, fill="#1f2933")
    # La columna de valores arranca tras la etiqueta mas larga del panel.
    col = x + 44 + max(ancho_texto(n, leyenda) for n, _, _ in leyendas) + 26
    for j, (nombre, color, valor) in enumerate(leyendas):
        yy = y + ch + 10 + j * 34
        d.rectangle((x + 10, yy + 8, x + 32, yy + 20), fill=color)
        d.text((x + 44, yy), nombre, font=leyenda, fill="#1f2933")
        d.text((col, yy), valor, font=leyenda, fill="#52606d")

d.text((22, H - 66), NOTA, font=pie, fill="#52606d", spacing=8)
img.save(FIGURAS + "/conformeros.png")
print(f"escrito {FIGURAS}/conformeros.png  {img.size}")
