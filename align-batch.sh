# ncbi_dataset contains NCBI reference Klebsiella genome
# number 4... contains Zanzibar blood sample sequences

# pipeline outlined at
# https://chatgpt.com/share/6706a40b-4580-8005-878a-27906b0c83c3
# although the alignment is totally done with RagTag or NumCmer (and the latter visualised)

REF_FILE="ncbi_dataset/ncbi_dataset/data/GCA_000240185.2/GCA_000240185.2_ASM24018v2_genomic.fna"

# Define the base directory and the filename pattern
BASE_DIR="./number4"
#BASE_DIR="./tmp"
PATTERN="*.fna" 

mkdir ragtag_out

# Loop through all files matching the pattern in the base directory and its subdirectories
find "$BASE_DIR" -type f -name "$PATTERN" | while read -r file; do
    # Extract the filename without the directory path
    filename=$(basename "$file")

    # Define the output file using the original filename (with .out extension, for example)
    sam_output_file="${filename}_aligned.sam"
    bam_output_file="${filename}_aligned.bam"
    bam_output_sorted_file="${filename}_aligned_sorted.bam"
    img_file="${filename}_alignment.png"

    # Perform your desired operation on the file (e.g., process it)
    echo "Processing $file..."

    # SAMs and BAMs don't play a role in the final version
    #minimap2 -x asm5 -t 4 -a "$REF_FILE" "$file" > "$sam_output_file"
    #samtools view -S -b "$sam_output_file" > "$bam_output_file"
    #samtools sort "$bam_output_file" -o "$bam_output_sorted_file"
    #samtools flagstat "$bam_output_sorted_file"

    # scaffold with RagTag
    ragtag.py scaffold "$file" "$REF_FILE"

    # just label RagTag output and pull into a consistent place
    # Loop through all files in the subdirectory
    for ragtagfile in "ragtag_output"/*; do
	# Check if it's a file (skip directories)
	if [ -f "$ragtagfile" ]; then
            # Extract the directory and the filename
            ragtagdirname=$(dirname "$ragtagfile")
            ragtagfilename=$(basename "$ragtagfile")

            # Create the new filename by appending the string
            new_filename="${filename}-${ragtagfilename}"

            # Rename the file
            mv "$ragtagfile" "ragtag_out/$new_filename"

            echo "Renamed $ragtagfile to ragtag_out/$new_filename"
	fi
    done

    # hard to visualise RagTag output? so visualise alignment from nucmer
    nucmer --prefix=alignment "$REF_FILE" "$file"
    mummerplot --png --layout --filter alignment.delta
    # automated script production fails; use manual gnuplot call instead
    gnuplot "out.gp"
    cp "out.png" "$img_file"

    echo "Output written to $output_file"
done

