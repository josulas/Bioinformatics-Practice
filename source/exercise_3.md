# EJERCICIO 3

Este ejercicio se encarga de alinear secuencias de ADN o proteínas utilizando el programa MUSCLE. Se toma un conjunto de secuencias como entrada y se genera un archivo de salida que contiene las secuencias alineadas, ya sea en formato FASTA o CLUSTAL.

## Descripción del Código

El código realiza las siguientes funciones:

1. **Lectura de Secuencias**: Carga secuencias de un archivo FASTA de origen y de un archivo FASTA de secuencias objetivo.
2. **Generación de un Archivo FASTA**: Escribe las secuencias leídas en un archivo FASTA temporal.
3. **Ejecución de MUSCLE**: Utiliza el programa MUSCLE para alinear las secuencias, generando dos archivos de salida:
   - Un archivo en formato FASTA (`aligned_sequences.fasta`)
   - Un archivo en formato CLUSTAL (`aligned_sequences.clustal`)
4. **Limpieza**: Elimina el archivo temporal utilizado para el alineamiento.

## Inputs

El script toma los siguientes parámetros de entrada desde la línea de comandos:

1. `path_to_source_fasta_file`: Ruta al archivo FASTA que contiene las secuencias de origen.
2. `path_to_target_sequences_fasta_file`: Ruta al archivo FASTA que contiene las secuencias objetivo que se desean alinear

## Outputs

  - Un archivo en formato FASTA (`aligned_sequences.fasta`)
  - Un archivo en formato CLUSTAL (`aligned_sequences.clustal`)