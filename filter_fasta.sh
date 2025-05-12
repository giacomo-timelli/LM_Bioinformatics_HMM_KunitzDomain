#!/usr/bin/bash

# Input FASTA file
input="pdb_kunitz_csv.fasta"
# Output FASTA file
output="pdb_kunitz_csv_filtered.fasta"

# Initialize variables
> "$output"
header=""
sequence=""

# Read the input file line by line
while IFS= read -r line; do
  if [[ $line == ">"* ]]; then
    # If a sequence was previously stored, evaluate its length
    if [[ -n $sequence ]]; then
      len=${#sequence}
      if (( len >= 45 && len <= 80 )); then
        echo "$header" >> "$output"
        echo "$sequence" >> "$output"
      fi
      sequence=""
    fi
    header="$line"  # Save the new header
  else
    sequence+="$line"  # Concatenate the sequence lines
  fi
done < "$input"

# Evaluate the final sequence at the end of the file
if [[ -n $sequence ]]; then
  len=${#sequence}
  if (( len >= 45 && len =< 80 )); then
    echo "$header" >> "$output"
    echo "$sequence" >> "$output"
  fi
fi
