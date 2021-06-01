#!/bin/bash
set -euo pipefail

STAG_MWC_OUTPUT_PATH=/home/groups/lu_kpmi/outputs

rm -f table.csv
rm -rf kreport_files
mkdir -p kreport_files

echo "sample_name,patient_id,treatment,time" > sample_metadata.csv

for f in $(ls $STAG_MWC_OUTPUT_PATH | grep LV) ; do # TODO: Enable processing of any sample, not only the ones in old format
    sample_name=$(echo $f | awk -F '_' '{print $1"_"$2"_"$3}')
    patient_id=$(echo $f | awk -F '_' '{print $1}')
    treatment=$(echo $f | awk -F '_' '{print $2}')
    time=$(echo $f | awk -F '_' '{print $3}')

    for s in $STAG_MWC_OUTPUT_PATH/$f/* ; do
        if [ -d $s/kraken2 ] ; then
			if [ -e $s/kraken2/1.kreport ] ; then
	            cp $s/kraken2/1.kreport kreport_files/$sample_name.kreport
			elif [ -e $s/kraken2/2.kreport ] ; then
				cp $s/kraken2/2.kreport kreport_files/$sample_name.kreport
			elif [ -e $s/kraken/$sample_name.kreport ] ; then
				cp $s/kraken2/$sample_name.kreport kreport_files/$sample_name.kreport
			fi
            break
        fi
    done

	echo $sample_name

    echo "$sample_name,$patient_id,$treatment,$time" >> sample_metadata.csv
done

echo "Kreports extracted"

BRACKEN_DB_PATH=$GROUP/databases/full_ref_bafp/database150mers.kmer_distrib

rm -rf bracken_1
rm -rf bracken_10
mkdir -p bracken_1
mkdir -p bracken_10

for f in ./kreport_files/* ; do
	name=$(basename ${f%.kreport})
	echo $name
	../../Bracken/src/est_abundance.py -i $f -k $BRACKEN_DB_PATH -o bracken_1/$name.bracken -t 1 -l G --out-report bracken_1/$name.bracken.kreport
	../../Bracken/src/est_abundance.py -i $f -k $BRACKEN_DB_PATH -o bracken_10/$name.bracken -t 10 -l G --out-report bracken_10/$name.bracken.kreport
done

echo "Bracken finished"

python3 ./calc_richness.py > richness.csv

echo "Richness calculated"

Rscript richness_summary.R

echo "Richness summary complete"

python3 ./pairwise.py

echo "Pairwise analysis done"
