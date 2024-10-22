# Ejercicio 4

Este ejercicio tiene como objetivo realizar un análisis de secuencias biológicas utilizando herramientas de EMBOSS. A partir de un archivo de entrada en formato GenBank, el script identifica las secuencias de lectura abiertas (ORFs), encuentra la secuencia más larga y luego realiza un análisis de dominio sobre dicha secuencia.

## Descripción del Código

El script utiliza las siguientes herramientas y funciones:

1. **getORFs**: Esta función llama a la herramienta `getorf` de EMBOSS para extraer las ORFs del archivo de entrada. El tamaño mínimo de las ORFs se puede ajustar mediante el parámetro `min_size_nucleotides`, cuyo valor por defecto es 1200 nucleótidos. Las secuencias resultantes se devuelven en formato FASTA.

2. **get_largest_seq**: Esta función toma una lista de secuencias (en formato FASTA) y determina cuál es la secuencia más larga.

3. **perform_domain_analysis**: Utiliza la herramienta `patmatmotifs` de EMBOSS para realizar un análisis de dominios sobre la secuencia más larga. Los resultados se almacenan en un archivo de salida especificado.

## Input

El script acepta un archivo de entrada en formato GenBank y puede recibir un argumento opcional para habilitar la salida detallada. El uso es el siguiente:

- `path_to_genbank_file`: Ruta al archivo GenBank que contiene las secuencias a analizar.
- `path_to_output_dir`: Ruta del directorio donde se guardarán los archivos de salida.
- `verbose` (opcional): Si se establece en `true`, proporciona salida detallada en la consola.

## Output

El resultado del script es un archivo llamado `motif_analysis.dbmotif`, que contiene los resultados del análisis de dominios. La ruta del archivo de salida se imprime en la consola al finalizar el script. Este será escrito en la carpeta `path_to_output_dir/`, y se se encuentra un archivo con el mismo nombre se reescribirá.
