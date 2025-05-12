#!/usr/bin/bash

# Input file from PDBefold
input="pdb_kunitz_clstr_rep.ali"

# Output file in Stockholm format
output="pdb_kunitz_clstr_rep_ali_str.sth"

# Add Stockholm header
echo "# STOCKHOLM 1.0" > "$output"

# Convert alignment to Stockholm without introducing empty lines
awk '
BEGIN {
    seqname = ""
    seq = ""
}
{
    gsub(/\r/, "")  # remove carriage return (if present)
    if (substr($1, 1, 1) == ">") {
        if (seqname != "") {
            print seqname, seq >> output
        }
        seqname = toupper(substr($1, 2))
        gsub(/^PDB:/, "", seqname)
        gsub(":", "_", seqname)
        seq = ""
    } else {
        seq = seq toupper($0)
    }
}
END {
    if (seqname != "") {
        print seqname, seq >> output
    }
}
' output="$output" "$input"

# Add end-of-file tag for Stockholm
echo "//" >> "$output"





# Build the HMM model from the aligned sequences
hmmbuild pdb_kunitz_clstr_rep_ali_str.hmm pdb_kunitz_clstr_rep_ali_str.sth


