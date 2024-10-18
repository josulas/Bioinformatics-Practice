# Ejercicio 1

Este ejercicio consiste en leer un archivo GenBank, extraer la traducción de la secuencia de nucleótidos, traducir todas las Regiones Abiertas de Lectura (ORFs) y guardar una de ellas en formato FASTA (elegimos aquella que se expresa en la practica).

## Descripción

El script está diseñado para realizar las siguientes operaciones:

1. **Cargar un archivo GenBank**: Se especifica el nombre del archivo GenBank a leer.
2. **Extraer la traducción**: Busca y extrae la secuencia de traducción de las características de tipo CDS en el archivo.
3. **Traducir ORFs**: Traduce todas las ORFs en la secuencia, tanto en la hebra directa como en la hebra inversa complementaria, considerando 6 marcos de lectura.
4. **Encontrar la ORF correspondiente**: Identifica la ORF que contiene la traducción extraída.
5. **Guardar la ORF en formato FASTA**: La ORF encontrada se guarda en un archivo FASTA, utilizando el ID de acceso como encabezado.

## Input

- **Archivo GenBank**: El archivo especificado debe ser otorgado y debe estar en formato GenBank.

## Output

- **Archivo FASTA**: Se guarda una ORF en formato FASTA con el ID de acceso en el encabezado. El nombre del archivo será `{accession_id}_orf{index}.fasta`.