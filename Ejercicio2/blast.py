"""
This script reads sequences from FASTA files, 
runs a local BLAST search for each sequence against a local database, 
and runs a remote BLAST search for each sequence against the SwissProt database.
The output of each BLAST search is saved in an XML file.
"""


import subprocess
import os
from Bio.Blast import NCBIWWW
from Bio import SeqIO


def read_sequences_from_files(sequence_files: list[str]) -> list[str]:
    """
    Reads all the sequences from the given files and returns them as a list of strings.
    """
    sequences = []
    for file_name in sequence_files:
        with open(file_name, 'r', encoding='UTF-8') as file:
            for seq_record in SeqIO.parse("infile.fasta", "fasta"):
                sequences.append(str(seq_record.seq))
    return sequences


def run_local_blast(sequences: list[str], db_path: str, output_file_name: str):
    """
    Run a local BLAST search for each sequence in the given list of sequences.
    """
    for i, sequence in enumerate(sequences):
        output = f"{output_file_name}_local_{i+1}.xml"
        with open(f"query_{i+1}.fasta", 'w', encoding='UTF-8') as temp_file:
            temp_file.write(sequence)
        blast_command = [
            "blastp",
            "-query", f"query_{i+1}.fasta",
            "-db", db_path,
            "-outfmt", "5",
            "-out", output
        ]
        subprocess.run(blast_command)
        print(f"Resultado BLAST local guardado en: {output}")


def run_remote_blast(sequences: list[str], output_file_name: str):
    """
    Run a remote BLAST search for each sequence in the given list of sequences.
    """
    for i, sequence in enumerate(sequences):
        output = f"{output_file_name}_remote_{i+1}.xml"
        result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence)
        with open(output, "w", encoding="UTF-8") as out_file:
            out_file.write(result_handle.read())
        result_handle.close()
        print(f"Resultado BLAST remoto guardado en: {output}")


def main():
    """
    Main function.
    """
    # Obtener todos los archivos de secuencia en el directorio actual
    sequence_files = [f for f in os.listdir() if f.endswith('.fasta')]
    # Leer las secuencias desde los archivos
    sequences = read_sequences_from_files(sequence_files)
    # Rutas y nombres
    db_path = "path_to_your_local_db/swissprot_db"
    output_file = "blast_output"
    # Ejecución BLAST local
    run_local_blast(sequences, db_path, output_file)
    # Ejecución BLAST remota
    run_remote_blast(sequences, output_file)

if __name__ == "__main__":
    main()
