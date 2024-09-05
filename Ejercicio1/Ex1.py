"""
Instructions:
- The transcription GenBank file is stored in the variable SEQ_FILENAME. Modiy this variable to
change which file is being read.
- Execute the script in the same directory as the GenBank file.
- The script will load the GenBank file, extract the translation (if present), 
translate all the ORFs and save them to a .txt file. 
It will then find the ORF that contains the extracted translation.
- The script will save the specified ORF in FASTA format with the accession ID in the header. 
To change which ORF is saved, modify the ORF_TO_SAVE variable.
"""

from Bio import SeqIO
from Bio.Seq import Seq

SEQ_FILENAME = "transcript-1.gb"
ORF_TO_SAVE = 2

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

    
    def load(self, filename: str):
        """
        Loads GenBank file and stores the nucleotide sequence and accession ID
        """
        try:
            self.filename = filename
            with open(filename, "r", encoding='UTF-8') as file:
                record = SeqIO.read(file, "genbank")
                self.sequence = str(record.seq)
                self.accession = record.id
                print(f'\nGenBank file "{filename}" was successfully loaded. Sequence length: {len(self.sequence)} bases.')
        except ValueError:
            print(f"\nError: '{filename}' is not a valid GenBank file.")
        except FileNotFoundError:
            print(f"\nError: File '{filename}' not found.")


    def pad_sequence(self, seq: str):
        """
        Pads the sequence with 'N' bases to ensure the length is divisible by 3
        """
        if len(seq) % 3 != 0:
            padding_length = 3 - (len(seq) % 3)
            seq += 'N' * padding_length
            print(f"\n{padding_length} 'N' bases added to sequence.")
        return seq


    def translate_orfs(self):
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


    def save_orfs_to_file(self, filename: str):
        """
        Save all 6 ORFs to a text file 
        Note: It was utilized as an intermediate check to compare the reading frames
        """
        orfs_list = self.translate_orfs()  # Ensure ORFs are translated
        with open(filename, "w", encoding='UTF-8') as file:
            for idx, orf in enumerate(orfs_list, start=1):
                file.write(f"ORF {idx}:\n{orf}\n\n")
        print(f"\nORFs were successfully saved to '{filename}'.")


    def extract_translation(self):
        """
        Extract the translation from the GenBank file as a reference
        """
        try:
            with open(self.filename, "r", encoding="UTF-8") as file:
                record = SeqIO.read(file, "genbank")
                for feature in record.features:
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        self.translation = feature.qualifiers["translation"][0]
                        print(f"\nTranslation sequence was successfully extracted from '{self.filename}'.")
                        break
        except FileNotFoundError:
            print(f"\nError: File '{self.filename}' not found.")
        except ValueError:
            print(f"\nError: '{self.filename}' is not a valid GenBank file.")


    def find_matching_orf(self):
        """
        Find the ORF that contains the extracted translation
        """
        if not self.translation:
            print("\nNo translation found in the GenBank file.")
            return None
        orfs_list = self.translate_orfs()  # Ensure ORFs are translated
        for idx, orf in enumerate(orfs_list, start=1):
            if self.translation in str(orf):
                print(f"\nMatching ORF found: ORF {idx} contains the translated sequence.")
                return idx
        print("\nNo matching ORF found.")
        return None


    def save_orf_as_fasta(self, orf_index: int):
        """
        Save the specified ORF in FASTA format with accession ID in the header
        """
        if not hasattr(self, 'accession'):
            print("\nNo accession ID found. Cannot create FASTA header.")
            return
        orfs_list = self.translate_orfs()
        if orf_index < 1 or orf_index > len(orfs_list):
            print(f"\nInvalid ORF index: {orf_index}. Must be between 1 and {len(orfs_list)}.")
            return
        orf_sequence = orfs_list[orf_index - 1]
        # Find the position of the first Methionine (start codon)
        start_index = str(orf_sequence).find('M')
        if start_index == -1:
            print(f"\nNo 'M' found in ORF{orf_index}. Writing full sequence.")
            start_index = 0
        else:
            print(f"\nStarting sequence from the first 'M' at position {start_index} in ORF{orf_index}.")

        trimmed_sequence = str(orf_sequence)[start_index:]

        if trimmed_sequence == self.translation:
            print("\nThe trimmed sequence matches the extracted translation.")
        else:
            print("\nThe trimmed sequence does NOT match the extracted translation.")

        filename = f"{self.accession}_orf{orf_index}.fasta"
        with open(filename, "w", encoding="UTF-8") as file:
            file.write(f">{self.accession} ORF{orf_index}\n{trimmed_sequence}\n")
        print(f"\nORF{orf_index} was successfully saved to '{filename}' starting from the first 'M'.")


genbank_seq = GenBankReader(SEQ_FILENAME)
genbank_seq.extract_translation()
matching_orf = genbank_seq.find_matching_orf()

if matching_orf:
    orfs = genbank_seq.translate_orfs()
    print(f"\nTranslated sequence in the matching ORF:\n{orfs[matching_orf - 1]}")
else:
    print("\nNo matching ORF was found.")

genbank_seq.save_orf_as_fasta(ORF_TO_SAVE)
