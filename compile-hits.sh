# this is, for the moment, specific to dataset 4
# lots of those sequences were in fact E coli. to pull the klebsiella set, go through the dnadiff alignment reports and extract a list of references for which UnalignedBases (with respect to Klebsiella reference) is less than 20%
grep "UnalignedBases" number4/outputs/*report | sed 's/[()%]/ /g' | awk '{if($3 < 20) print $0;}' | awk 'BEGIN{FS=".";}{print $1;}' | awk 'BEGIN{FS="/";}{print $(NF);}' > kp-hits.txt

# now we want to pairwise align these klebsiella hits, and produce a report of ANI scores

# Define the directory where your reference genomes are stored
REFERENCE_DIR="number4/"
EXTENSION="fna"  # Change this if needed

# Read the list of references
references=()
while IFS= read -r line; do
    references+=("$line")
done < kp-hits.txt

# initialise a summary file
echo "" > summary.txt

# Loop over all pairs of references
for ((i=0; i<${#references[@]}; i++)); do
    for ((j=i+1; j<${#references[@]}; j++)); do
        ref1="${references[i]}"
        ref2="${references[j]}"

        # Run dnadiff
        echo "Comparing $ref1 vs $ref2..."
        dnadiff  \
            "${REFERENCE_DIR}/${ref1}_annotation/${ref1}.${EXTENSION}" \
            "${REFERENCE_DIR}/${ref2}_annotation/${ref2}.${EXTENSION}"
	cp out.report "${ref1}-${ref2}.report"
	# extract the first (1-to-1) AvgIdentity score (ANI) and append to a summary
	grep "AvgIdentity" "${ref1}-${ref2}.report" | awk -v v1="$ref1" -v v2="$ref2" 'BEGIN{n=0;}{if(n==0) { print v1, v2, $3; } n++;}' >> summary.txt

    done
done
