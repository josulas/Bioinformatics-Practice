"""
Instructions can be found in 'exercise_1.md'
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq


class GenBankReader():
    """
    Class to read a GenBank file and perform various operations on the sequence data.
    """

    def __init__(self, filename: str|None = None):
        # Initialize class
        self.sequence = ""
        self.translation = ""
        self.accession = ""
        self.orfs = []
        if filename:
            self.load(filename)

    
    def load(self, filename: str, verbose: bool = False):
        """
        Loads GenBank file and stores the nucleotide sequence and accession ID
        """
        try:
            self.filename = filename
            with open(filename, "r", encoding='UTF-8') as file:
                record = SeqIO.read(file, "genbank")
                self.sequence = str(record.seq)
                self.accession = record.id
                if verbose:
                    print(f'\nGenBank file "{filename}" was successfully loaded. Sequence length: {len(self.sequence)} bases.')
        except ValueError:
                print(f"\nError: '{filename}' is not a valid GenBank file.")
        except FileNotFoundError:
                print(f"\nError: File '{filename}' not found.")


    def pad_sequence(self, seq: str, verbose: bool = False):
        """
        Pads the sequence with 'N' bases to ensure the length is divisible by 3
        """
        if len(seq) % 3 != 0:
            padding_length = 3 - (len(seq) % 3)
            seq += 'N' * padding_length
            if verbose:
                print(f"\n{padding_length} 'N' bases added to sequence.")
        return seq


    def translate_orfs(self, verbose: bool = False):
        """
        Translate nucleotide sequence in all 3x2=6 ORFs (forward and reverse)
        """
        sequence = Seq(self.sequence)

        if not self.orfs:  # Only translate if ORFs are not already computed
            translated_orfs = []

            # Forward strand
            for i in range(3):
                frame = sequence[i:]
                padded_frame = self.pad_sequence(frame)
                orf = Seq(padded_frame).translate(to_stop=True)
                translated_orfs.append(orf)
                if verbose:
                    print(f"\nForward strand reading frame {i+1} was successfully translated.")

            # Reverse complement strand
            reverse_complement = sequence.reverse_complement()
            for i in range(3):
                frame = reverse_complement[i:]
                padded_frame = self.pad_sequence(frame)
                orf = Seq(padded_frame).translate(to_stop=True)
                translated_orfs.append(orf)
                if verbose:
                    print(f"\nReverse strand reading frame {i+1} was successfully translated.")

            self.orfs = translated_orfs  # Store the translated ORFs
        return self.orfs


    def save_orfs_to_file(self, filename: str, verbose: bool = False):
        """
        Save all 6 ORFs to a text file 
        Note: It was utilized as an intermediate check to compare the reading frames
        """
        orfs_list = self.translate_orfs()  # Ensure ORFs are translated
        with open(filename, "w", encoding='UTF-8') as file:
            for idx, orf in enumerate(orfs_list, start=1):
                file.write(f"ORF {idx}:\n{orf}\n\n")
        if verbose:
            print(f"\nORFs were successfully saved to '{filename}'.")


    def extract_translation(self, verbose = False):
        """
        Extract the translation from the GenBank file as a reference
        """
        try:
            with open(self.filename, "r", encoding="UTF-8") as file:
                record = SeqIO.read(file, "genbank")
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        self.translation = feature.qualifiers["translation"][0]
                        if verbose: 
                            print(f"\nTranslation sequence was successfully extracted from '{self.filename}'.")
                        break
        except FileNotFoundError:
            print(f"\nError: File '{self.filename}' not found.")
        except ValueError:
            print(f"\nError: '{self.filename}' is not a valid GenBank file.")


    def find_matching_orf(self, verbose: bool = False):
        """
        Find the ORF that contains the extracted translation
        """
        if not self.translation:
            print("\nNo translation found in the GenBank file.")
            return None
        orfs_list = self.translate_orfs()  # Ensure ORFs are translated
        for idx, orf in enumerate(orfs_list, start=1):
            if self.translation in str(orf):
                if verbose: 
                    print(f"\nMatching ORF found: ORF {idx} contains the translated sequence.")
                return idx
        if verbose: 
            print("\nNo matching ORF found.")
        return None


    def save_orf_as_fasta(self, orf_index: int, dir: str = r"./",verbose: bool = False):
        """
        Save the specified ORF in FASTA format with accession ID in the header
        """
        if not hasattr(self, 'accession'):
            print("\nNo accession ID found. Cannot create FASTA header.") if verbose else None
            return
        orfs_list = self.translate_orfs()
        if orf_index < 1 or orf_index > len(orfs_list):
            print(f"\nInvalid ORF index: {orf_index}. Must be between 1 and {len(orfs_list)}.") if verbose else None
            return
        orf_sequence = orfs_list[orf_index - 1]
        # Find the position of the first Methionine (start codon)
        start_index = str(orf_sequence).find('M')
        if start_index == -1:
            print(f"\nNo 'M' found in ORF{orf_index}. Writing full sequence.") if verbose else None
            start_index = 0
        else:
            print(f"\nStarting sequence from the first 'M' at position {start_index} in ORF{orf_index}.") if verbose else None

        trimmed_sequence = str(orf_sequence)[start_index:]

        if trimmed_sequence == self.translation:
            print("\nThe trimmed sequence matches the extracted translation.") if verbose else None
        else:
            print("\nThe trimmed sequence does NOT match the extracted translation.") if verbose else None

        filename = os.path.join(dir, f"{self.accession}_orf{orf_index}.fasta")
        with open(filename, "w", encoding="UTF-8") as file:
            file.write(f">{self.accession} ORF{orf_index}\n{trimmed_sequence}\n")
        print(f"\nORF{orf_index} was successfully saved to '{filename}' starting from the first 'M'.") if verbose else None
        return filename


if __name__ == "__main__":
    error_msg = "Usage: python3 exercise_1.py <path_to_genbank_file> <path_to_output_dir> (optional) <verbose: true | false>"
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print(error_msg)
        sys.exit(1)
    if len(sys.argv) == 4:
        verbose_order = sys.argv[3]
        if verbose_order == 'true':
            verbose = True
        elif verbose_order == 'false':
            verbose = False
        else:
            print(error_msg)
            sys.exit(1)
    else:
        verbose = False
    genbank_file_name = sys.argv[1]
    output_dir = sys.argv[2]
    genbank_seq = GenBankReader(genbank_file_name)
    genbank_seq.extract_translation(verbose=verbose)
    matching_orf = genbank_seq.find_matching_orf(verbose=verbose)
    if matching_orf:
        orfs = genbank_seq.translate_orfs(verbose=verbose)
        if verbose:
            print(f"\nTranslated sequence in the matching ORF:\n{orfs[matching_orf - 1]}")
        output_file = genbank_seq.save_orf_as_fasta(matching_orf, output_dir, verbose=verbose)
        print(output_file)
        sys.exit(0)
    else:
        if verbose:
            print("\nNo matching ORF was found.")
        sys.exit(1)
