from Bio.Blast import NCBIWWW, NCBIXML
import os

def read_sequences_from_files(sequence_files):
    sequences = []
    for file_name in sequence_files:
        with open(file_name, 'r') as file:
            sequences.append(file.read())
    return sequences

def run_local_blast(sequences, db_path, output_file):
    for i, sequence in enumerate(sequences):
        output = f"{output_file}_local_{i+1}.xml"
        with open(f"query_{i+1}.fasta", 'w') as temp_file:
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

def run_remote_blast(sequences, output_file):
    for i, sequence in enumerate(sequences):
        output = f"{output_file}_remote_{i+1}.xml"
        result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence)
        with open(output, "w") as out_file:
            out_file.write(result_handle.read())
        result_handle.close()
        print(f"Resultado BLAST remoto guardado en: {output}")

def main():
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
