
# Ejercicio 2

Este script realiza búsquedas BLAST locales y remotas para secuencias de ADN o proteínas extraídas de archivos FASTA. Se ejecutan las siguientes operaciones:

1. **Lectura de Secuencias**: Se leen secuencias desde un archivo en formato FASTA.

2. **Búsqueda BLAST Local**:
   - Se ejecuta una búsqueda BLAST local para cada secuencia contra una base de datos local (previamente descargada).
   - El resultado de la búsqueda se guarda en un archivo XML. Si el archivo ya existe, no se realiza la búsqueda.

3. **Búsqueda BLAST Remota**:
   - Se ejecuta una búsqueda BLAST remota contra la base de datos SwissProt utilizando los servidores de NCBI.
   - El resultado de esta búsqueda también se guarda en un archivo XML. Si el archivo ya existe, no se realiza la búsqueda.

4. **Extracción de Hits**:
   - A partir de los resultados en XML, se extraen las 10 mejores secuencias encontradas (basado en el puntaje).
   - Estas secuencias se guardan en archivos en formato FASTA.

## Requisitos

- Es necesario un sistema operativo **Linux** con y **Biopython** y **BLAST** instalados.

## Input

1. `path_to_fasta_file`: El nombre de un archivo FASTA que contiene una o más secuencias de ADN o proteínas.
2. `path_to_output_dir`: Ruta del directorio donde se guardarán los archivos de salida.

## Output

- Archivos **XML** que contienen los resultados de las búsquedas BLAST (tanto local como remota).
- Archivos **FASTA** con las mejores secuencias (hits) obtenidas en ambas búsquedas (local y remota).
- Los archivos serán escritos en la carpeta `path_to_output_dir/`, y se se encuentra un archivo con el mismo nombre, este se reescribirá.
