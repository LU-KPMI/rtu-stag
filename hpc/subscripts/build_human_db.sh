#!/bin/bash
#PBS -N build_db
#PBS -l feature=largescratch
#PBS -l nodes=1:ppn=32,pmem=6g
#PBS -l walltime=96:00:00
#PBS -q long
#PBS -j oe
#PBS -d .

module load conda
source /opt/exp_soft/conda/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate stag-mwc

work_dir=/scratch/kc01/human_genome_build

rm -rf $work_dir
mkdir -p $work_dir

cp -r ../../kraken2 $work_dir
cp human_genome.fasta $work_dir

cd $work_dir

kraken2/kraken2-build --download-taxonomy --db human_genome --use-ftp
kraken2/kraken2-build --add-to-library human_genome.fasta --db human_genome
kraken2/kraken2-build --build --db human_genome --threads 32
