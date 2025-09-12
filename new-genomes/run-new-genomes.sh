# updated 12 sep 2025
chmod +x *.sh

# unpack new genome FASTA files
cd number138
tar -xzf new-genomes.tar.gz
cd ..

## do alignments of number138 directory with Kp reference
./align-batch.sh 138
## extract genuine Kp hits, fold in with old Olav data, and build script comparing pairwise everywhere 
./process-pairs.sh 138

chmod +x part_*
./part_00 > tmp00 &
./part_01 > tmp01 &
./part_02 > tmp02 &
./part_03 > tmp03 &
./part_04 > tmp04 &
./part_05 > tmp05 &
./part_06 > tmp06 &
./part_07 > tmp07 &
./part_08 > tmp08 &
./part_09 > tmp09 &
./part_10 > tmp10 &
