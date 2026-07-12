"""
Este script calcula los contactos de las estructuras PDB filtradas utilizando la herramienta PLIP.
Se asume que los archivos PDB filtrados se encuentran en la carpeta "structures/filtered" 
y los resultados se guardarán en la carpeta "structures/contacts" solo manteniendo los archivos XML.
"""

import os
import subprocess   

INPUT_FOLDER = "structures/filtered"
OUTPUT_FOLDER = "structures/contacts"

skips = [
    "7AAA",
    "9DMC"
]


# Ejecutar la función de descarga si el script se ejecuta directamente
if __name__ == "__main__":
    # Crear la carpeta de salida si no existe
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Directorio creado: {OUTPUT_FOLDER}")

    pdb_files = os.listdir(INPUT_FOLDER)
    print(f"Se han encontrado {len(pdb_files)} estructuras para analizar:\n")

    for code in pdb_files:

        if code[:4] in skips:
            print(f"✗ Se omite {code[:4]} (en la lista de skips)")
            continue

        input_path = os.path.join(INPUT_FOLDER, code)
        output_path = os.path.join(OUTPUT_FOLDER, code[:4])

        if os.path.exists(output_path):
            print(f"  -> {code[:4]} ya procesado. Se omite.")
            continue

        os.makedirs(output_path, exist_ok=True)

        print(f"Calculando contactos para {code[:4]}...")

        # Llamada al script externo para calcular contactos
        
        cmd = [
                "plip",
                "-f", str(input_path),
                "-o", str(output_path),
                "-x",
                "-t"
            ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print(f"   ✓ Contactos calculados correctamente para {code[:4]}, borrando archivos intermedios...")

            for archivo in os.listdir(output_path):
                if not archivo.endswith(".xml"):
                    os.remove(os.path.join(output_path, archivo))
                
        else:
            print(f"   ✗ Error al calcular contactos para {code[:4]}")
            print(result.stderr)