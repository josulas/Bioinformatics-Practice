#!/bin/bash

# Create log directory, if it does not exists
log_directory="./log/"
if [ ! -d $log_directory ]; then
    mkdir $log_directory
fi

# Create current log file
log_file=$log_directory$(date +"%Y-%m-%d_%T").log
touch $log_file

# Create the output directory
output_directory="./output/"
if [ ! -d $output_directory ]; then
    mkdir $output_directory
fi

# run makeconfig.sh
./source/makeconfig.sh
exit_status_config=$?
if [ $exit_status_config -eq 0 ]; then
    echo "Configuration completed successfully" >> $log_file
else
    echo "Configuration failed with exit code $exit_status_config" >> $log_file
    exit $exit_status_config
fi

# Exercise 1
# List all .gb files in the input directory
input_directory='input/'
gb_files=($(ls $input_directory*.gb 2>/dev/null))
# Check if there are any .gb files
if [ ${#gb_files[@]} -eq 0 ]; then
    echo "No GenBank (.gb) files found in the input directory"
    echo "No input files (.gb) were found" >> $log_file
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
    echo "Invalid choice. Exiting"
    echo "User made an invalid choice on the selection of the .gb file for exercise 1" >> $log_file
    exit 1
fi
# Get the selected genbank file
selected_file_1="${gb_files[$((choice-1))]}"
echo "You selected: $selected_file_1"
# Call excercise_1.py to find the ORFs
output_file_1=$(python3.10 source/exercise_1.py "$selected_file_1" "$output_directory")
exit_status_1=$?
if [ $exit_status_1 -eq 0 ]; then
    echo "Python script for exercise 1 completed successfully" >> $log_file
    echo "Output file: $output_file_1" >> $log_file
    echo "The output for the first exercise has been saved at $output_file_1"
else
    echo "Python script for exercise 1 failed with exit code $exit_status_1" >> $log_file
    echo "An error occurred during the execution of the first exercise. Exiting"
    exit $exit_status_1
fi

# Exercise 2
output_files_2_string=$(python3.10 source/exercise_2.py "$output_file_1" "$output_directory")
IFS=' ' read -r -a output_files_2 <<< "$output_files_2_string"
exit_status_2=$?
if [ $exit_status_2 -eq 0 ]; then
    echo "Python script for exercise 2 completed successfully" >> $log_file
    echo "Output files: $output_files_2_string" >> $log_file
    echo "The output for the second exercise has been saved at $output_files_2_string"
else
    echo "Python script for exercise 2 failed with exit code $exit_status_2" >> $log_file
    echo "An error occurred during the execution of the second exercise. Exiting"
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
    echo "User made an invalid choice on the selection of the .fasta file for the exercise 3" >> $log_file
    exit 1
fi
selected_file_3="${output_files_2[$((choice_3 - 1))]}"
output_files_3_string=$(python3.10 source/exercise_3.py "$output_file_1" "$selected_file_3" "$output_directory")
exit_status_3=$?
if [ $exit_status_3 -eq 0 ]; then
    echo "Python script for exercise 3 completed successfully" >> $log_file
    echo "Output files: $output_files_3_string" >> $log_file
    echo "The MSA results for the third exercise has been saved at $output_files_3_string"
else
    echo "Python script for exercise 3 failed with exit code $exit_status_3" >> $log_file
    echo "An error occurred during the execution of the third exercise. Exiting"
    exit $exit_status_3
fi

# Exercise 4
output_file_4=$(python3.10 source/exercise_4.py "$selected_file_1" "$output_directory")
exit_status_4=$?
if [ $exit_status_4 -eq 0 ]; then
    echo "Python script for exercise 4 completed successfully" >> $log_file
    echo "Output file: $output_file_4" >> $log_file
    echo "The report for the fourth exercise has been saved at $output_file_4"
else
    echo "Python script for exercise 4 failed with exit code $exit_status_4" >> $log_file
    echo "An error occurred during the execution of the fourth exercise. Exiting"
    exit $exit_status_4
fi

#Exercise 5
output_file_5=$(python3.10 source/exercise_5.py "$selected_file_1" "$output_directory")
exit_status_5=$?
if [ $exit_status_5 -eq 0 ]; then
    echo "Python script for exercise 5 completed successfully" >> $log_file
    echo "Output files: $output_file_5" >> $log_file
    echo "The list of primers for the fifth exercise has been saved at $output_file_5"
else
    echo "Python script for exercise 5 failed with exit code $exit_status_5" >> $log_file
    echo "An error occurred during the execution of the fifth exercise. Exiting"
    exit $exit_status_5
fi