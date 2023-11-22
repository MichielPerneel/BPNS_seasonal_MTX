#!/bin/bash

# Run this script as create_phyloDB_extended_annotation_file.sh FASTA_FILE TEXT_FILE
# where FASTA_FILE is the path to the input fasta file and TEXT_FILE is the path to the output text file.

# Check if the number of arguments is correct
if [ $# -ne 2 ]; then
  echo "Usage: $0 FASTA_FILE TEXT_FILE"
  exit 1
fi

# Store the fasta file and text file as variables
fasta_file=$1
text_file=$2

# Extract the headers from the fasta file
headers=$(grep "^>" $fasta_file)

# Initialize an empty text file
> $text_file

# Prompt the user for the text to be added to the additional columns
read -p "Enter text for column 2: " col2
read -p "Enter text for column 3: " col3
read -p "Enter text for column 4: " col4

# Loop through the headers, adding them to the text file
# with the text from the additional columns
while read -r line; do
    modified_header=${line:1}
    echo -e "$modified_header\t$col2\t$col3\t$col4" >> $text_file
done <<< "$headers"
