#!/bin/bash

work_path="$PWD/../../"
taxon_db_path="/home/groups/lu_kpmi/databases/full_ref_bafp"
human_ref_path="/home/groups/lu_kpmi/databases/human_reference"
sample_path="/home/groups/lu_kpmi/renamed_samples"
resistome_path="/home/groups/lu_kpmi/databases/groot_db/arg-annot_index"
output_path="/home/groups/lu_kpmi/outputs"

for f in ${sample_path}/*_1.fq.gz; do # for name generation, don't want to trigger twice - limiting myself to the first file of the pair
    name=$(echo $f | sed 's:.*/::' | sed 's/_[^_]*$//')
    read_1=$sample_path/${name}_1.fq.gz
    read_2=$sample_path/${name}_2.fq.gz
    qsub subscripts/sub.run.sh -F "$name $read_1 $read_2 $work_path $taxon_db_path $human_ref_path $resistome_path $output_path" # create jobs for all of the samples
done
