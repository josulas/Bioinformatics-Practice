# Blast configuration
# Define the expected database files
db_folder="blast_database"
db_files=("swissprot" "swissprot_db.phr" "swissprot_db.pin" "swissprot_db.psq")
# Function to check if all database files are present
check_db_files() {
    if [ ! -d "$db_folder" ]; then
        return 1
    fi
    cd $db_folder
    for file in "${db_files[@]}"; do
        if [ ! -f "$file" ]; then
            cd ..
            return 1  # Return 1 if any file is missing
        fi
    done
    cd ..
    return 0  # Return 0 if all files are present
}
# Check if the swissprot database is already downloaded and formatted
if (! check_db_files ); then
    echo "SwissProt database not found. Downloading and formatting..."
    if [ ! -d "$db_folder" ]; then
        mkdir $db_folder
    fi
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
    gunzip swissprot.gz
    /bin/makeblastdb -in swissprot -dbtype prot -out swissprot_db
    mv swissprot* $db_folder
    # Recheck if the files were successfully created
    if check_db_files; then
        echo "SwissProt database downloaded and formatted successfully."
    else
        echo "Error: Failed to download and format the SwissProt database."
        echo "Exiting"
        exit 1
    fi
fi

# Prosite Configuration
prosite_directory="prosite"
if [ ! -d "$prosite_directory" ]; then
    mkdir $prosite_directory
    wget https://ftp.expasy.org/databases/prosite/prosite.doc -P $prosite_directory --quiet
    if [ ! $? -eq 0 ]; then
        echo "wget couldn't fetch prosite.doc file. Exiting"
        exit 1
    fi
    wget https://ftp.expasy.org/databases/prosite/prosite.dat -P $prosite_directory --quiet
    if [ ! $? -eq 0 ]; then
        echo "wget couldn't fetch prosite.dat file. Exiting" 
        exit 1
    fi
    echo admin | sudo -S prosextract -prositedir $prosite_directory -verbose N
    if [ ! $? -eq 0 ]; then
        echo "prosextract couldn't set up prosite. Exiting" 
        exit 1
    fi
fi