#!/bin/bash
set -euo pipefail

STAG_MWC_OUTPUT_PATH=/home/groups/lu_kpmi/outputs_reformat

rm -f table.csv
rm -rf bracken_files

mkdir -p bracken_files

echo "sample_name,patient_id,treatment,time" > table.csv

for f in $(ls $STAG_MWC_OUTPUT_PATH | grep LV) ; do
    sample_name=$(echo $f | awk -F '_' '{print $1"_"$2"_"$3}')
    patient_id=$(echo $f | awk -F '_' '{print $1}')
    treatment=$(echo $f | awk -F '_' '{print $2}')
    time=$(echo $f | awk -F '_' '{print $3}')

    for s in $STAG_MWC_OUTPUT_PATH/$f/* ; do
        if [ -d $s/kraken2 ] ; then
            cp $s/kraken2/*_bracken_species.kreport bracken_files/$sample_name.kreport
            break
        fi
    done

    echo "$sample_name,$patient_id,$treatment,$time" >> table.csv
done

echo "Bracken outputs extracted"

python3 ../../kraken-biom/kraken_biom.py bracken_files/*.kreport -o abundance.biom --fmt json

echo "Created biom files"

Rscript richness.R

echo "Richness calculated"
