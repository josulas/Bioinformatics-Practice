
# Ejercicio 2

Este script realiza búsquedas BLAST locales y remotas para secuencias de ADN o proteínas extraídas de archivos FASTA. Se ejecutan las siguientes operaciones:

1. **Lectura de Secuencias**: Se leen secuencias desde un archivo en formato FASTA.

2. **Búsqueda BLAST Local**:
   - Se ejecuta una búsqueda BLAST local para cada secuencia contra una base de datos local (previamente descargada).
   - El resultado de la búsqueda se guarda en un archivo XML.

3. **Búsqueda BLAST Remota**:
   - Se ejecuta una búsqueda BLAST remota contra la base de datos SwissProt utilizando los servidores de NCBI.
   - El resultado de esta búsqueda también se guarda en un archivo XML.

4. **Extracción de Hits**:
   - A partir de los resultados en XML, se extraen las 10 mejores secuencias encontradas (basado en el puntaje).
   - Estas secuencias se guardan en archivos en formato FASTA.

## Requisitos

- El script debe ejecutarse desde la misma carpeta donde están presentes los archivos FASTA.
- Es necesario un sistema operativo **Linux** con y **BLAST** instalados.

## Input

- Un archivo FASTA que contiene una o más secuencias de ADN o proteínas.

## Output

- Archivos **XML** que contienen los resultados de las búsquedas BLAST (tanto local como remota).
- Archivos **FASTA** con las mejores secuencias (hits) obtenidas en ambas búsquedas (local y remota).
