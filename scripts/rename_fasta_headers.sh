#!/bin/bash

# Small script to rename fasta headers
# Run this script like:
# ./rename_fasta_headers.sh PREFIX input.fasta > output.fasta

# Check if the user provided a string to use as the prefix for the new headers
if [ -z "$1" ]
then
    echo "Error: No prefix provided."
    exit 1
fi

# Check if the user provided a fasta file
if [ -z "$2" ]
then
    echo "Error: No fasta file provided."
    exit 1
fi

# Initialize the contig count to zero
count=0

# Read the fasta file and rename the headers
while read -r line
do
    # If the line starts with ">" it is a header line
    if [[ $line =~ ^\> ]]
    then
        # Increment the contig count
        count=$((count+1))

        # Print the new header line
        echo ">$1_$count"
    else
        # If it is not a header line, just print it as is
        echo "$line"
    fi
done < "$2"