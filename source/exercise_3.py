"""
Instructions can be found in 'exercise_3.md'
"""

import sys
import subprocess
import os
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def write_sequences_to_fasta(sequences: List[str], ids: List[str], fasta_file):
    """
    Writes a list of sequences to a FASTA file.
    
    Args:
        sequences (list of tuples): A list of (sequence_id, sequence_string) tuples.
        fasta_file (str): The output FASTA filename.
    """
    records = []
    for i, seq in enumerate(sequences):
        record = SeqRecord(Seq(seq), id=ids[i], description="")
        records.append(record)
    
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


def run_muscle(input_fasta: str, output_file: str, output_mode: str = 'fasta', verbose: bool = False):
    if output_mode == 'fasta':
        muscle_cline = [
            "muscle",                   # MUSCLE executable
            "-in", input_fasta,         # Input FASTA
            "-out", output_file,        # Output FASTA with aligned sequences
            '-quiet'
        ]
    elif output_mode == 'clustal':
        muscle_cline = [
            "muscle",                   # MUSCLE executable
            "-in", input_fasta,         # Input FASTA
            "-out", output_file,        # Output CLUSTAL with aligned sequences
            '-clw', '-quiet'
        ]
    else:
        raise ValueError(f"output_mode must be 'fasta' or 'clustal', got {output_mode} instead.")
    try:
        subprocess.run(muscle_cline, check=True)
        print(f"Alignment saved to: {output_file}") if verbose else None
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE: {e}")
        sys.exit(1)


if __name__ == "__main__":
    error_message = "Usage: python3 exercise_3.py <path_to_source_fasta_file> <path_to_target_sequences_fasta_file> <path_to_output_dir> (optional) <verbose: true | false>"
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print()
        sys.exit(1)
    if len(sys.argv) == 5:
        verbose_order = sys.argv[4]
        if verbose_order == 'true':
            verbose = True
        elif verbose_order == 'false':
            verbose = False
        else:
            print(error_message)
            sys.exit(1)
    else:
        verbose = False
    source_seqs_file = sys.argv[1]
    target_seqs_file = sys.argv[2]
    output_dir = sys.argv[3]

    process_id = lambda seq: "".join(seq.description.split(";")[0].replace("RecName: Full=", "").replace(" ", ""))
    sequence_list = []
    id_list = []

    # Get which sequence should be chosen from the source_seqs_file
    digits = []
    for char in target_seqs_file:
        if char.isdigit():
            digits.append(char)
    source_seq_indx = int("".join(digits)) - 1

    with open(source_seqs_file, "r") as handle:
        source_seq = list(SeqIO.parse(handle, "fasta"))[source_seq_indx]
    sequence_list.append(str(source_seq.seq))
    id_list.append(process_id(source_seq))

    with open(target_seqs_file, "r") as handle:
        target_seqs = list(SeqIO.parse(handle, "fasta"))
    for seq in target_seqs:
        sequence_list.append(str(seq.seq))
        id_list.append(process_id(seq))

    input_fasta = ".msa_input"
    output_fasta = os.path.join(output_dir, "aligned_sequences.fasta")
    output_clustal = os.path.join(output_dir, "aligned_sequences.clustal")
    output_files = [output_fasta, output_clustal]

    write_sequences_to_fasta(sequence_list, id_list, input_fasta)
    # Run MUSCLE to align the sequences
    run_muscle(input_fasta, output_fasta)
    run_muscle(input_fasta, output_clustal, 'clustal')
    os.remove(input_fasta)
    print(" ".join(output_files))
