#!/bin/bash
set -euo pipefail

STAG_MWC_OUTPUT_PATH=/home/groups/lu_kpmi/outputs

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

rm -rf bracken
mkdir -p bracken

for f in ./kreport_files/* ; do
	name=$(basename ${f%.kreport})
	echo $name
	../../Bracken/src/est_abundance.py -i $f -k $BRACKEN_DB_PATH -o bracken/$name.bracken -t 1 -l G --out-report bracken/$name.bracken.kreport
done

echo "Bracken finished"

python3 ./gen_table.py > abundances.csv

echo "Abundance table generated"

python3 ./calc_richness.py > richness.csv

echo "Richness calculated"

python3 ./graphs.py

echo "Richness summary visualization complete"

#python3 ./pairwise.py
#
#echo "Pairwise analysis done"
