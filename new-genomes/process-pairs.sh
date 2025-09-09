# Define the base directory and the filename pattern
if [ "$1" == "1" ]; then
  BASE_DIR="./number1"
fi
if [ "$1" == "4" ]; then
  BASE_DIR="./number4"
fi
if [ "$1" == "138" ]; then
  BASE_DIR="./number138"
fi

PATTERN="*.fasta"
# Define the directory where your reference genomes are stored
REFERENCE_DIR=$BASE_DIR
EXTENSION="fasta"  # Change this if needed

# lots of those sequences were in fact E coli. to pull the klebsiella set, go through the dnadiff alignment reports and extract a list of references for which UnalignedBases (with respect to Klebsiella reference) is less than 20%
grep "UnalignedBases" ${BASE_DIR}/outputs/*report | sed 's/[()%]/ /g' | awk '{if($3 < 20) print $0;}'  | awk 'BEGIN{FS="fasta";}{print $1;}' | awk 'BEGIN{FS="/";}{print $(NF);}' | sed 's/[.]//g' | sort | uniq > $1-kp-hits.txt

# create scripts for an all-against-all comparison combining new Kp data and old

# get lists of references and filenames
ls From_Olav_fixed/*.fasta > old-kp-files.txt
sed 's/From_Olav_fixed[/]//g' old-kp-files.txt | sed 's/[.]fasta//g' > old-kp.txt
awk '{printf("./number138/%s.fasta\n", $0);}' 138-kp-hits.txt > 138-kp-hits-files.txt
cat 138-kp-hits.txt old-kp.txt > all-kp.txt
cat 138-kp-hits-files.txt old-kp-files.txt > all-kp-files.txt

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
	    'BEGIN{printf("dnadiff -p %s-%s %s %s\n", ref1, ref2, fname1, fname2);
printf("rm %s-%s.1coords\nrm %s-%s.1delta\nrm %s-%s.delta\nrm %s-%s.mcoords\n", ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2);
printf("rm %s-%s.mdelta\nrm %s-%s.qdiff\nrm %s-%s.rdiff\nrm %s-%s.snps\nrm %s-%s.unqry\nrm %s-%s.unref\nrm %s-%s.mgaps\nrm %s-%s.ntref\n", ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2, ref1, ref2);}' >> compare-script.sh
    done
done

# split this big script into many parts for parallelisation
csplit -f part_ -n 2 compare-script.sh 37000 {9}
