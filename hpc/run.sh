#!/bin/bash

datestamp=$(date -d "today" +"%Y%m%d%H%M")

mkdir -p $HOME/cron_logs
LOG_FILE=$HOME/cron_logs/$datestamp

sample_list="list" # The pathname of file containing list of samples to be processed. Each line should be in following format: sample_name path_to_read_1 path_to_read_2
work_path="" # The pathname of working directory you created
taxon_db_path="/home/groups/lu_kpmi/databases/full_ref_bafp"
human_ref_path="/home/groups/lu_kpmi/databases/human_reference"
resistome_path="/home/groups/lu_kpmi/databases/groot_db/arg-annot_index"
output_path="/home/groups/lu_kpmi/outputs"

echo "Contents of file:" >> $LOG_FILE
cat $sample_list >> $LOG_FILE

echo "Start processing..." >> $LOG_FILE
touch $sample_list.tmp

ENABLE_QC_READS=True
ENABLE_HOST_REMOVAL=True
ENABLE_KRAKEN2=True
ENABLE_GROOT=True
ENABLE_AMRPLUSPLUS=True
BRACKEN_TRESH=1

while read sample; do
    name=$(echo $sample | awk -F ' ' '{print $1}')
    read_1=$(echo $sample | awk -F ' ' '{print $2}')
    read_2=$(echo $sample | awk -F ' ' '{print $3}')

    qsub $work_path/rtu-stag/hpc/subscripts/sub.run.sh -F "$name $read_1 $read_2 $work_path $taxon_db_path $human_ref_path $resistome_path $output_path \
$ENABLE_QC_READS $ENABLE_HOST_REMOVAL $ENABLE_KRAKEN2 $ENABLE_GROOT $ENABLE_AMRPLUSPLUS $BRACKEN_TRESH"
    if [ $? -eq 0 ]; then
        echo "$name put in queue" >> $LOG_FILE
    else
        echo "$name NOT put in queue" >> $LOG_FILE
        echo $sample >> $sample_list.tmp
    fi
done < $sample_list

mv $sample_list.tmp $sample_list
