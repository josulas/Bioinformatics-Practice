#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3

class GenBankReader:

    def __init__(self, filename=None):
        # Initialize class
        self.sequence = ""
        self.translation = ""
        self.accession = ""
        self.orfs = []  
        if filename:
            self.load(filename)

    def load(self, filename):
        # Loads GenBank file and stores the nucleotide sequence and accession ID
        try:
            self.filename = filename
            with open(filename, "r") as file:
                record = SeqIO.read(file, "genbank")
                self.sequence = str(record.seq)
                self.accession = record.id 
                print(f'\nGenBank file "{filename}" was successfully loaded. Sequence length: {len(self.sequence)} bases.')
        except Exception as e:
            print(f"\nError: {e}")

    def pad_sequence(self, seq):
        # Add trailing 'N's to ensure the sequence length is divisible by 3
        if len(seq) % 3 != 0:
            padding_length = 3 - (len(seq) % 3)
            seq += 'N' * padding_length
            print(f"\n{padding_length} 'N' bases added to sequence.")
        return seq

    def translate_orfs(self):
        # Translate nucleotide sequence in all 3x2 ORFs (forward and reverse)
        sequence = Seq(self.sequence)

        if not self.orfs:  # Only translate if ORFs are not already computed
            translated_orfs = []

            # Forward strand
            for i in range(3):
                frame = sequence[i:]
                padded_frame = self.pad_sequence(frame)
                orf = Seq(padded_frame).translate(to_stop=True)
                translated_orfs.append(orf)
                print(f"\nForward strand reading frame {i+1} was successfully translated.")

            # Reverse complement strand
            reverse_complement = sequence.reverse_complement()
            for i in range(3):
                frame = reverse_complement[i:]
                padded_frame = self.pad_sequence(frame)
                orf = Seq(padded_frame).translate(to_stop=True)
                translated_orfs.append(orf)
                print(f"\nReverse strand reading frame {i+1} was successfully translated.")

            self.orfs = translated_orfs  # Store the translated ORFs
        return self.orfs

    def save_orfs_to_file(self, filename):
        # Save all 6 ORFs to a text file - was utilized as an intermediate check to compare the reading frames
        try:
            orfs = self.translate_orfs()  # Ensure ORFs are translated
            with open(filename, "w") as file:
                for idx, orf in enumerate(orfs, start=1):
                    file.write(f"ORF {idx}:\n{orf}\n\n")
            print(f"\nORFs were successfully saved to '{filename}'.")
        except Exception as e:
            print(f"\nError: {e}")

    def extract_translation(self):
        # Extract the translation from the GenBank file as a reference
        try:
            with open(self.filename, "r") as file:
                record = SeqIO.read(file, "genbank")
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        self.translation = feature.qualifiers["translation"][0]
                        print(f"\nTranslation sequence was successfully extracted from '{self.filename}'.")
                        break
        except Exception as e:
            print(f"\nError: {e}")

    def find_matching_orf(self):
        # Find the ORF that contains the extracted translation
        if not self.translation:
            print("\nNo translation found in the GenBank file.")
            return None

        orfs = self.translate_orfs()  # Ensure ORFs are translated
        for idx, orf in enumerate(orfs, start=1):
            if self.translation in str(orf):
                print(f"\nMatching ORF found: ORF {idx} contains the translated sequence.")
                return idx
        print("\nNo matching ORF found.")
        return None

    def save_orf_as_fasta(self, savename, orf_index):
        # Save the specified ORF in FASTA format with accession ID in the header
        try:
            orfs = self.translate_orfs()
            if orf_index < 1 or orf_index > len(orfs):
                print(f"\nInvalid ORF index: {orf_index}. Must be between 1 and {len(orfs)}.")
                return
            
            orf_sequence = orfs[orf_index - 1]

            # Find the position of the first Methionine (start codon)
            start_index = str(orf_sequence).find('M')
            if start_index == -1:
                print(f"\nNo 'M' found in ORF{orf_index}. Writing full sequence.")
                start_index = 0
            else:
                print(f"\nStarting sequence from the first 'M' at position {start_index} in ORF{orf_index}.")

            trimmed_sequence = str(orf_sequence)[start_index:]

            if not hasattr(self, 'accession'):
                print("\nNo accession ID found. Cannot create FASTA header.")
                return

            if trimmed_sequence == self.translation:
                print("\nThe trimmed sequence matches the extracted translation.")
            else:
                print("\nThe trimmed sequence does NOT match the extracted translation.")

            filename = f"{self.accession}_orf{orf_index}.fasta"
            with open(filename, "w") as file:
                file.write(f">{self.accession} ORF{orf_index}\n{trimmed_sequence}\n")
            print(f"\nORF{orf_index} was successfully saved to '{filename}' starting from the first 'M'.")


        except Exception as e:
            print(f"\nError: {e}")            


# Using the script 
genbank_seq = GenBankReader('transcript-1.gb')
genbank_seq.extract_translation()
matching_orf = genbank_seq.find_matching_orf()

if matching_orf:
    orfs = genbank_seq.translate_orfs()
    print(f"\nTranslated sequence in the matching ORF:\n{orfs[matching_orf - 1]}")
else:
    print("\nNo matching ORF was found.")

genbank_seq.save_orf_as_fasta('orf2.fasta', orf_index=2)
