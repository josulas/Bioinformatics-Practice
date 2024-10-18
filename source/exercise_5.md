# Ejercicio 5

Este script genera primers a partir de una secuencia de ADN en formato GenBank, evaluando diferentes parámetros como el contenido de GC y la temperatura de melting (Tm). El script también asegura que los primers no comiencen ni terminen con las bases G o C. Los parámetros para la generación de los primers se leen desde un archivo de configuración en formato JSON.


## Descripción

1. **Lectura de secuencias**: El script lee una secuencia desde un archivo en formato GenBank usando la función `read_sequences_from_file`. Solo se utiliza la primera secuencia encontrada.

2. **Cálculo del contenido de GC**: A través de la función `gc_content`, el script calcula el porcentaje de guanina (G) y citosina (C) en una secuencia de ADN.

3. **Cálculo de la temperatura de melting (Tm)**: Utilizando la fórmula de Wallace, la función `melting_temperature` calcula la temperatura de melting de la secuencia.

4. **Validación de extremos**: La función `valid_ends` verifica que el primer no comience ni termine con las bases G o C.

5. **Generación de primers**: La función `generate_primers` genera posibles primers a partir de la secuencia de ADN, considerando las restricciones establecidas en el archivo de configuración. Los parámetros de configuración incluyen:
   - Longitud mínima y máxima del primer.
   - Porcentaje mínimo y máximo de contenido de GC.
   - Temperatura máxima de melting (Tm).
   
6. **Archivo de configuración**: El archivo `parameters_5.json` debe contener los siguientes parámetros:
   ```json
   {
       "primer_length_min": <longitud mínima>,
       "primer_length_max": <longitud máxima>,
       "gc_min": <porcentaje mínimo de GC>,
       "gc_max": <porcentaje máximo de GC>,
       "melting_temp_max": <temperatura máxima de melting>
   }

## Input

- **Archivo GenBank**: El script espera como primer argumento la ruta a un archivo en formato GenBank que contiene la secuencia de ADN

- **Archivo de configuración JSON**: El script lee un archivo llamado parameters_5.json que contiene los parámetros necesarios para la generación de primers, como longitud mínima y máxima, contenido de GC y temperatura de melting.

## Output

- **Archivo txt**: El script genera un archivo de texto llamado primers.txt que contiene una lista de los primers generados.