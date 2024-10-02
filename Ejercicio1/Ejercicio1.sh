#!/bin/bash

# Habilitar modo estricto para manejar errores
set -e


# Mensaje de inicio
echo "Ejecutando el Ejercicio 1: procesamiento de archivo GenBank."

# Verificar si Python está instalado
if ! command -v python &> /dev/null; then
    echo "Error: Python no está instalado. Por favor, instálalo."
    exit 1
fi

# Crear un entorno virtual si no existe
if [ ! -d "venv" ]; then
    echo "Creando entorno virtual..."
    python -m venv venv
fi

# Activar el entorno virtual
source venv/Scripts/activate  # Cambia esto para Windows

# Instalar dependencias si es necesario
if [ ! -f "requirements.txt" ]; then
    echo "Creando requirements.txt..."
    echo "biopython" > requirements.txt
fi

pip install -r requirements.txt

# Ejecutar el script de Python con las variables
echo "Ejecutando el script de Python para analizar el archivo  y guardar ORF."

python Ex1.py 

# Confirmar la finalización
echo "Proceso completado."