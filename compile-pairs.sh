# create scripts for an all-against-all comparison combining new Kp data and old

# get lists of references and filenames
ls From_Olav_fixed/*.fasta > old-kp-files.txt
sed 's/From_Olav_fixed[/]//g' old-kp-files.txt | sed 's/[.]fasta//g' > old-kp.txt
awk '{printf("./number1/%s.fasta\n", $0);}' 1-kp-hits.txt > 1-kp-hits-files.txt
cat 1-kp-hits.txt old-kp.txt > all-kp.txt
cat 1-kp-hits-files.txt old-kp-files.txt > all-kp-files.txt

# Read the list of references
references=()
while IFS= read -r line; do
    references+=("$line")
done < all-kp.txt

# Read the list of filenames
fnames=()
while IFS= read -r line; do
    fnames+=("$line")
done < all-kp-files.txt

# initialise a summary file
echo "" > "compare-script.sh"

# Loop over all pairs of references
for ((i=0; i<${#references[@]}; i++)); do
    for ((j=i+1; j<${#references[@]}; j++)); do
        ref1="${references[i]}"
        ref2="${references[j]}"
	fname1="${fnames[i]}"
	fname2="${fnames[j]}"

	awk -v ref1="$ref1" -v ref2="$ref2" -v fname1="$fname1" -v fname2="$fname2" \
	    'BEGIN{printf("dnadiff -p %s-%s %s %s\n", ref1, ref2, fname1, fname2);}' >> compare-script.sh
    done
done

# split this big script into many parts for parallelisation
csplit -f part_ -n 2 compare-script.sh 2000 {9}
