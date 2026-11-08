# Compone cada fotograma: imagen de ChimeraX + gráfico de RMSD + título, reloj y aviso de cámara lenta.
import csv
import os
import sys

from PIL import Image, ImageDraw, ImageFont

# Por defecto, la raiz del repositorio.
ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
D = ROOT + "/run/render/"           # en .gitignore
TABLAS = ROOT + "/results/md/tables"
os.makedirs(D + "final", exist_ok=True)
FUENTE = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
NEGRITA = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
titulo = ImageFont.truetype(NEGRITA, 30)
reloj = ImageFont.truetype(FUENTE, 26)
aviso = ImageFont.truetype(NEGRITA, 22)

for p in csv.DictReader(open(TABLAS + "/plan_render_60fps.csv")):
    n = int(p["salida"])
    img = Image.open(f"{D}mol/{n:04d}.png").convert("RGBA")
    inset = Image.open(f"{D}inset/{n:04d}.png").convert("RGBA")
    x0, y0 = 32, img.height - inset.height - 20
    panel = Image.new("RGBA", img.size, (0, 0, 0, 0))
    ImageDraw.Draw(panel).rounded_rectangle((x0 - 8, y0 - 8, x0 + inset.width + 8, y0 + inset.height + 8),
                                            radius=14, fill=(255, 255, 255, 235), outline="#d9dee5", width=2)
    img.alpha_composite(panel)
    img.alpha_composite(inset, (x0, y0))
    d = ImageDraw.Draw(img)
    d.text((40, 32), "PARP1 · ligando 453 · 50 ns de dinámica molecular", font=titulo, fill="#1f2933")
    d.text((40, 76), f"t = {float(p['tiempo_ns']):.1f} ns", font=reloj, fill="#52606d")
    if p["lento"] == "1":
        d.rounded_rectangle((36, 118, 600, 156), radius=8, fill="#fde8d7")
        d.text((48, 123), "cámara lenta ×6 · cambio de modo de unión", font=aviso, fill="#c2410c")
    img.convert("RGB").save(f"{D}final/{n:04d}.png")
