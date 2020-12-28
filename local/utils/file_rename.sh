#!/bin/bash

# future develop.
# show error message instructing the user on correct usage and create an examplary CSV file.
# rsync the files straight to HPC already renamed and, if required, put together.
pkgs='pv'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
  # sudo apt install $pkgs
  exit "pipe-viewer not installed, please run 'sudo apt install pv' to continue"
fi

mkdir -p renamed
for i in $(awk -F "," '{OFS=","; print $10, $11}' sample_list_data.csv); do # Read csv file.
  if [[ "${i}" == *\-* ]]; # check if sample spans multiple barcodes (has a dash in old name). If has a dash perform seq with comma delim.
  then
    a=$(echo ${i} | awk -F "_" '{sub(/-.*$/, "", $3); print $3 }') # get barcode range start.
    b=$(echo ${i} | awk -F "," '{sub(/.*-/, "", $1); print $1 }') # get barcode range stop.
    barcodes=$(seq -s , $a $b) 
    oldname=$(echo ${i} | awk -F "," '{print $1 }') # get old name from csv.
    newname=$(echo ${i} | awk -F "," '{print $2 }') # get new name from csv.
    # merge and rename R1
    # 'eval' not the best solution can introduce known vulnerability, but is allowed in this case, should rewrite in future.
    searchname=${oldname%_$a-*}
    eval ls "${searchname}_{${barcodes}}_1.fq.gz" | xargs pv | cat > "renamed/${newname}_1.fq.gz" | echo "${newname}_1.fq.gz" # rename R1 file combining seq and csv variables.
    # merge and rename R2
    eval ls "${searchname}_{${barcodes}}_2.fq.gz" | xargs pv | cat > "renamed/${newname}_2.fq.gz" | echo "${newname}_2.fq.gz"  # rename R2 file combining seq and csv variables.
  else # rename the files that dont span multiple barcodes by copying to avoid potential catastrophe.
    oldname2=$(echo ${i} | awk -F "," '{print $1 }') # get old name from csv.
    newname2=$(echo ${i} | awk -F "," '{print $2 }') # get new name from csv. 
    # rename R1
    rsync --progress "${oldname2}_1.fq.gz" "renamed/${newname2}_1.fq.gz" 
    # rename R2
    rsync --progress "${oldname2}_2.fq.gz" "renamed/${newname2}_2.fq.gz" 
  fi
done
find renamed/*fq.gz -size 0 -delete # a workaround to get rid of empty files. Should avoid making empty files in the first place.