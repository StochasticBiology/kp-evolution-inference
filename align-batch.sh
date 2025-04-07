#!/bin/bash

# ncbi_dataset contains NCBI reference Klebsiella genome
# number 4... contains Zanzibar blood sample sequences

# pipeline outlined at
# https://chatgpt.com/share/6706a40b-4580-8005-878a-27906b0c83c3
# although the alignment is totally done with RagTag or NumCmer (and the latter visualised)

#REF_FILE="ncbi_dataset/ncbi_dataset/data/GCA_000240185.2/GCA_000240185.2_ASM24018v2_genomic.fna"
REF_FILE="GCA_000240185.2_ASM24018v2_genomic.fna"

# Define the base directory and the filename pattern
if [ "$1" == "1" ]; then
  BASE_DIR="./number1"
  PATTERN="*.fasta"
fi
if [ "$1" == "4" ]; then
  BASE_DIR="./number4"
  PATTERN="*.fasta"
fi
if [ "$1" == "138" ]; then
  BASE_DIR="./number138"
  PATTERN="*.fasta"
fi

# allocate and create output directory
OUT_DIR="${BASE_DIR}/outputs" 
mkdir "$OUT_DIR"

# Loop through all files matching the pattern in the base directory and its subdirectories
find "$BASE_DIR" -type f -name "$PATTERN" | while read -r file; do
    # Extract the filename without the directory path
    filename=$(basename "$file")

    report_file="${filename}.report"

    # Perform your desired operation on the file (e.g., process it)
    echo "Processing $file..."
    # get report of alignment with respect to reference sequence
    echo "Running dnadiff"
    dnadiff "$REF_FILE" "$file"
    cp "out.report" "${OUT_DIR}/${report_file}"
    
    echo "Output written to $output_file"
done

