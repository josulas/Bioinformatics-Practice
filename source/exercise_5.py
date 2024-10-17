import sys
from typing import List
from Bio import SeqIO


def read_sequences_from_file(sequence_file: str) -> List[str]:
    """
    Reads all the sequences from the given files and returns them as a list of strings.
    """
    sequences = []
    with open(sequence_file, 'r', encoding='UTF-8') as file:
        for seq_record in SeqIO.parse(file, "genbank"):
            sequences.append(str(seq_record.seq))
    return sequences


def gc_content(sequence):
    """Calcula el porcentaje de GC en una secuencia."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def melting_temperature(sequence):
    """Calcula la temperatura de melting (Tm) usando la fórmula de Wallace."""
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    return 4 * gc_count + 2 * at_count


def valid_ends(sequence):
    """Verifica que la secuencia no tenga G o C en los extremos."""
    return not (sequence.startswith('G') or sequence.startswith('C') or 
                sequence.endswith('G') or sequence.endswith('C'))


def generate_primers(transcript_sequence, config):
    """Genera primers a partir de la secuencia del transcripto utilizando los parámetros del archivo de configuración."""
    primers = []
    length_min = config['primer_length_min']
    length_max = config['primer_length_max']
    gc_min = config['gc_min']
    gc_max = config['gc_max']
    melting_temp_max = config['melting_temp_max']
    
    for start in range(0, len(transcript_sequence) - length_min + 1):
        for length in range(length_min, length_max + 1):
            if start + length > len(transcript_sequence):
                break
            primer_candidate = transcript_sequence[start:start + length]
            gc_percentage = gc_content(primer_candidate)
            tm = melting_temperature(primer_candidate)
            
            if (gc_min <= gc_percentage <= gc_max and
                tm <= melting_temp_max and
                valid_ends(primer_candidate)):
                primers.append({
                    'sequence': primer_candidate,
                    'gc_content': gc_percentage,
                    'tm': tm
                })
                
    return primers

import json

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 exercise_5.py <path_to_fasta_file>")
        sys.exit(1)
    sequence_file_name = sys.argv[1]
    # Leer las secuencias desde el archivo
    transcript_sequence = read_sequences_from_file(sequence_file_name)[0]

    # Leer archivo de configuración
    with open('parameters_5.json') as f:
        config = json.load(f)
    
    # Generar primers
    primers = generate_primers(transcript_sequence, config)
    
    # Guardar resultado
    output_file_name = "primers.txt"
    with open(output_file_name, 'w') as file:
        for primer in primers:
            file.write(f"Primer: {primer['sequence']}, GC: {primer['gc_content']:.2f}%, Tm: {primer['tm']}°C\n")
    print(output_file_name)

if __name__ == "__main__":
    main()