"""
This script reads sequences from FASTA files, 
runs a local BLAST search for each sequence against a local database, 
and runs a remote BLAST search for each sequence against the SwissProt database.
The output of each BLAST search is saved in an XML file. Pay attention to:
- The script should be run from its folder
- .fasta files should be present in the folder in order to generate an output
- A linux machine is required, with biopython and blast installed
"""

import subprocess
import os
import sys
from typing import List
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def parse_blast_xml(xml_file: str, output_fasta: str, top_n: int = 10, verbose: bool = False):
    hits = []
    with open(xml_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Capture the score, title, and sequence for each hit
                    hit = {
                        'score': hsp.score,
                        'title': alignment.title,
                        'sequence': hsp.sbjct  # subject sequence
                    }
                    hits.append(hit)
    hits.sort(key=lambda x: x['score'], reverse=True)
    top_hits = hits[:top_n]
    records = []
    for i, hit in enumerate(top_hits):
        seq_record = SeqRecord(Seq(hit['sequence']), id=f"hit_{i+1}", description=hit['title'])
        records.append(seq_record)
    SeqIO.write(records, output_fasta, "fasta")
    if verbose:
        print(f"Top {top_n} sequences written to {output_fasta}")


def read_sequences_from_file(sequence_file: str) -> List[str]:
    """
    Reads all the sequences from the given files and returns them as a list of strings.
    """
    sequences = []
    with open(sequence_file, 'r', encoding='UTF-8') as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            sequences.append(str(seq_record.seq))
    return sequences


def run_local_blast(sequence: str, db_path: str, output_file_name: str, verbose: bool = False):
    """
    Run a local BLAST search for each sequence in the given list of sequences.
    """
    temp_file_name = ".fasta"
    with open(temp_file_name, 'w', encoding='UTF-8') as temp_file:
        temp_file.write(sequence)
    blast_command = [
        "blastp",
        "-query", temp_file_name,
        "-db", db_path,
        "-outfmt", "5",
        "-out", output_file_name
    ]
    subprocess.run(blast_command)
    if verbose:
        print(f"Resultado BLAST local guardado en: {output_file_name}")
    os.remove(temp_file_name)


def run_remote_blast(sequence: str, output_file_name: str, verbose: bool = False):
    """
    Run a remote BLAST search for each sequence in the given list of sequences.
    """
    result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence)
    with open(output_file_name, "w", encoding="UTF-8") as out_file:
        out_file.write(result_handle.read())
    result_handle.close()
    if verbose:
        print(f"Resultado BLAST remoto guardado en: {output_file_name}")


if __name__ == "__main__":
    # Rutas y nombres
    db_path = "./blast_database/swissprot_db"
    output_file_common = "blast_output"
    best_seq_fasta_file = "bests"
    # Obtenemos el nombre del archivo desde sys
    if len(sys.argv) < 2:
        print("Usage: python3 exercise_2.py <path_to_fasta_file>")
        sys.exit(1)
    sequence_file_name = sys.argv[1]
    # Leer las secuencias desde el archivo
    sequences = read_sequences_from_file(sequence_file_name)
    # Ejecutamos blast y guardamos el resultado
    return_files = []
    for index, sequence in enumerate(sequences):
        # Ejecución BLAST local
        output_file_local = f"{output_file_common}_local_{index + 1}.xml"
        if not os.path.exists(output_file_local):
            run_local_blast(sequence, db_path, output_file_local)
        # Guardamos los top 10 hits de la búsqueda local
        fasta_file_local = f"{best_seq_fasta_file}_local_{index + 1}.fasta"
        parse_blast_xml(output_file_local, fasta_file_local)
        return_files.append(fasta_file_local)
        # Ejecución BLAST remota
        output_file_remote = f"{output_file_common}_remote_{index + 1}.xml"
        if not os.path.exists(output_file_remote):
            run_remote_blast(sequence, output_file_remote)
        # Guardamos los top 10 hits de la búsqueda remota
        fasta_file_remote = f"{best_seq_fasta_file}_remote_{index + 1}.fasta"
        parse_blast_xml(output_file_remote, fasta_file_remote)
        return_files.append(fasta_file_remote)
    print(" ".join(return_files), end='')
    