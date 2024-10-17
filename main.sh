#!/bin/bash

# Exercise 1
# List all .gb files in the current directory
gb_files=($(ls *.gb 2>/dev/null))
# Check if there are any .gb files
if [ ${#gb_files[@]} -eq 0 ]; then
    echo "No GenBank (.gb) files found in the current directory."
    exit 1
fi
# Display the numbered menu
echo "Please select a GenBank file to process:"
for i in "${!gb_files[@]}"; do
    printf "%s) %s\n" "$((i+1))" "${gb_files[$i]}"
done
# Read user choice
read -p "Enter the number of the file: " choice_1
# Validate input
if [[ "$choice_1" -lt 1 || "$choice_1" -gt "${#gb_files[@]}" ]]; then
    echo "Invalid choice. Exiting."
    exit 1
fi
# Get the selected genbank file
selected_file_1="${gb_files[$((choice-1))]}"
echo "You selected: $selected_file_1"
# Call excercise_1.py to find the ORFs
output_file_1=$(python3.10 source/exercise_1.py "$selected_file_1")
exit_status_1=$?
if [ $exit_status_1 -eq 0 ]; then
    echo "Python script for exercise 1 completed successfully."
    echo "Generated file: $output_file_1"
else
    echo "Python script for exercise 1 failed with exit code $exit_status_1"
    echo "Exiting"
    exit $exit_status_1
fi

# Exercise 2
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
    mkdir $db_folder
    ./source/load_and_prepare_db.sh
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
output_files_2_string=$(python3.10 source/exercise_2.py "$output_file_1")
IFS=' ' read -r -a output_files_2 <<< "$output_files_2_string"
exit_status_2=$?
if [ $exit_status_2 -eq 0 ]; then
    echo "Python script for exercise 2 completed successfully."
    echo "Best hits saved at: $output_files_2_string"
else
    echo "Python script for exercise 2 failed with exit code $exit_status_2"
    echo "Exiting"
    exit $exit_status_2
fi

# Exercise 3
# Display the numbered menu
echo "Please select a .fasta file to align with a sequence in ${output_file_1}:"
for i in "${!output_files_2[@]}"; do
    printf "%s) %s\n" "$((i+1))" "${output_files_2[$i]}"
done
# Read user choice
read -p "Enter the number of the file: " choice_3
# Validate input
if [[ "$choice_3" -lt 1 || "$choice_3" -gt "${#output_files_2[@]}" ]]; then
    echo "Invalid choice. Exiting."
    exit 1
fi
selected_file_3="${output_files_2[$((choice_3 - 1))]}"
python3.10 source/exercise_3.py "$output_file_1" "$selected_file_3"
exit_status_3=$?
if [ $exit_status_3 -eq 0 ]; then
    echo "Python script for exercise 3 completed successfully."
    # echo "Generated file: $output_file_1"
else
    echo "Python script for exercise 3 failed with exit code $exit_status_3"
    echo "Exiting"
    exit $exit_status_3
fi

# Exercise 4
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
output_file_4=$(python3.10 source/exercise_4.py "$selected_file_1")
exit_status_4=$?
if [ $exit_status_4 -eq 0 ]; then
    echo "Python script for exercise 4 completed successfully."
    echo "Report created in: $output_file_4"
else
    echo "Python script for exercise 4 failed with exit code $exit_status_4"
    echo "Exiting"
    exit $exit_status_4
fi

#Exercise 5
output_file_5=$(python3.10 source/exercise_5.py "$selected_file_1")
exit_status_5=$?
if [ $exit_status_5 -eq 0 ]; then
    echo "Python script for exercise 5 completed successfully."
    echo "Primers saved in: $output_file_5"
else
    echo "Python script for exercise 5 failed with exit code $exit_status_5"
    echo "Exiting"
    exit $exit_status_5
fi