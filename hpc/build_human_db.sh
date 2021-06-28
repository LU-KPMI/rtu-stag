#!/bin/bash

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=568336000&extrafeat=null&conwithfeat=on&hide-cdd=on" -O chry.fasta # From https://www.ncbi.nlm.nih.gov/nuccore/CM000686.2?report=fasta
sed -i '1c>chrY' chry.fasta

wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz"
gunzip "chm13.draft_v1.1.fasta.gz"

cat chm13.draft_v1.1.fasta chry.fasta > human_genome.fasta

sed -i '/>/ s/$/|kraken:taxid|9606/' human_genome.fasta

rm "chm13.draft_v1.1.fasta" "chry.fasta"

qsub subscripts/build_human_db.sh
