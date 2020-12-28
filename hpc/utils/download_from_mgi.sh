#!/bin/bash
#PBS -N dl_mgi
#PBS -l nodes=1:ppn=1,pmem=6g
#PBS -l walltime=96:00:00
#PBS -q long
#PBS -j oe

# DO NOT USE, NOT COMPLETE!
# Logging into MGI server
lftp ftp://LU_metagenome@10.245.1.138
# A workaround
echo "set ssl:verify-certificate no"
# Mirror folder straight to raw_mgi_data folder, should check if this mirrors, without deleting from destination
lftp LU_metagenome@10.245.1.138:/> mirror -c --use-pget-n=10  --parallel=2  home/ /mnt/home/groups/lu_kpmi/raw_mgi_data


