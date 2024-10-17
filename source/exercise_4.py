import sys
import subprocess
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq


PROSITE_DAT_FILE_PATH = "prosite.dat"


def getORFs(input_file_name: str, min_size_nucleotides: int = 1200):
    emboss_command = [
        "getorf",
        "-minsize", str(min_size_nucleotides),
        "-sequence", input_file_name,
        "-outseq", "stdout"
    ]
    completion = subprocess.run(emboss_command, capture_output=True, text=True)
    if completion.returncode:
        print("Emboss getorf call was unsuccessfull")
        exit(1)
    sequences = SeqIO.parse(StringIO(completion.stdout), "fasta")
    return sequences


def get_largest_seq(sequences) -> str:
    largest = ""
    for sequence in sequences:
        if len(largest) < len(sequence.seq):
            largest = sequence.seq
    return str(largest)


def perform_domain_analysis(input_sequence: str, output_file_name: str):
    emboss_command = [
        "patmatmotifs",   # The EMBOSS tool for motif search
        "-sequence", "stdin",  # Provide sequence via stdin
        "-outfile", output_file_name,  # Specify where to store the results
        "-full",  # Include full details of the motifs
    ]
    completion = subprocess.run(emboss_command, input=input_sequence, capture_output=True, text=True)
    if completion.returncode:
        print("Emboss patmatmotifs call was unsuccessfull")
        exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 exercise_4.py <path_to_genbank_file> (optional) <verbose: true | false>")
        sys.exit(1)
    if len(sys.argv) >= 3:
        order = sys.argv[2]
        if order == 'true':
            verbose = True
        elif order == 'false':
            verbose = False
        else:
            print("Usage: python3 exercise_4.py <path_to_genbank_file> (optional) <verbose: true | false>")
            sys.exit(1)
    else:
        verbose = False
    output_file_path = "motif_analysis.dbmotif"
    genbank_file_name = sys.argv[1]
    seqs = getORFs(genbank_file_name)
    most_likely_peptide = get_largest_seq(seqs)
    perform_domain_analysis(most_likely_peptide, output_file_path)
    print(output_file_path)