# Read the list of references
references=()
while IFS= read -r line; do
    references+=("$line")
done < all-kp.txt

echo "" > summary.txt

# Loop over all pairs of references
for ((i=0; i<${#references[@]}; i++)); do
    for ((j=i+1; j<${#references[@]}; j++)); do
        ref1="${references[i]}"
        ref2="${references[j]}"
	fname1="${fnames[i]}"
	fname2="${fnames[j]}"

	grep "AvgIdentity" "${ref1}-${ref2}.report" | awk -v v1="$ref1" -v v2="$ref2" 'BEGIN{n=0;}{if(n==0) { print v1, v2, $3; } n++;}' >> summary.txt
    done
done
