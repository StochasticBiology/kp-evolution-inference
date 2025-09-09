#!/bin/bash
cd ../clean

for file in *.tsv; do
    echo Doing $file
    if [[ $(wc -l "$file" | cut -f1 -d" ") -eq 2 ]]
    then
        echo Only one sample, skipping
        continue
    fi
    python3 ../trees/LINcoding.py -i "$file" -t 1 -l 2- -b 3.02,7,69.79,93.16,98.41,98.89,99.36,99.68,99.84,100 -tree Y -lg Y
    mv ARBRE_REF.txt "${file%.*}".nwk
    echo done
done

cd -
