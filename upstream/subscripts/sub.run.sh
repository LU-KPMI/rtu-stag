#!/bin/bash
#PBS -N run_sample
#PBS -l feature=largescratch
#PBS -l nodes=1:ppn=32,pmem=6g
#PBS -l walltime=24:00:00
#PBS -q long
#PBS -j oe

# how many threads do we have?
threads=32
run_humann=false # NB - this takes a while to run

module load conda
# a bit of a stupid solution - but if it works it works
source /opt/exp_soft/conda/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate stag-mwc

name="$1" # Both reads will be renamed to $name_1.fq.gz and $name_2.fq.gz
read_1="$2" # Pathname for 1st read
read_2="$3" # Pathname for 2nd read
work_path="$4"
taxon_db_path="$5"
human_ref_path="$6"
resistome_path="$7"
output_path="$8"
export ENABLE_QC_READS=$9
export ENABLE_HOST_REMOVAL=${10}
export ENABLE_KRAKEN2=${11}
export ENABLE_GROOT=${12}
export ENABLE_AMRPLUSPLUS=${13}
export BRACKEN_TRESH=${14}
export prefix="/scratch/$(whoami)"
f="${prefix}/$(whoami)_$name"

echo "On machine $HOSTNAME"
echo "Start time" $(date)
echo "Sample $name"

# Use scratch dir to keep us from nuking their network infrastructure
rm -rf "$f" # clear out the folder in case this sample has already been on this nod
mkdir -p "$f"

{
    flock 200 # Multiple jobs running on the same node can start copying database, therefore add a lock
    # # copy the database folders over - just use scratch instead of using the sample dir
    # rm -rf "${prefix}/databases"
    if [ ! -d "${prefix}/databases" ]; then # NB: this will cause issues if we ever want to update the databases
        mkdir "${prefix}/databases"
        cp -r "${human_ref_path}" "${prefix}/databases"
        cp -r "${taxon_db_path}" "${prefix}/databases"
        cp -r "${resistome_path}" "${prefix}/databases"
    fi
} 200>$work_path/db_lock

# Copy stag, samples and config file
cp -r "${work_path}/stag-mwc" "$f"
envsubst < "${work_path}/rtu-stag/upstream/config.yaml" > "$f/stag-mwc/config.yaml"
mkdir "$f/stag-mwc/input"
cp $read_1 "$f/stag-mwc/input/${name}_1.fq.gz"
cp $read_2 "$f/stag-mwc/input/${name}_2.fq.gz"

# Copy kraken2
cp -r "${work_path}/kraken2" "${prefix}"

# Launch stag
cd "$f/stag-mwc"
flock $prefix/conda_lock snakemake --use-conda --create-envs-only --conda-prefix=$prefix/conda_envs
snakemake --use-conda --conda-prefix=$prefix/conda_envs --cores $threads || exit 1 # Exit if stag fails

cd ${prefix}
if [ "$run_humann" = true ] ; then
    # run the humann2 stuff outside of stag - just ripping the whole thing to deal with dep conflicts between humann2 and snakemake
    humann2_dir="$f/stag-mwc/output_dir/humann2/"
    metaphlan_dir="$f/stag-mwc/output_dir/metaphlan/"
    mkdir -p "$humann2_dir"
    mkdir -p "$metaphlan_dir"
    # at this point we know that host_removal samples exist due to them being made for groot
    echo "#SampleID\t$sample" > mpa2_table-v2.7.7.txt
    # metaphlan had to run before humann2
    metaphlan --input_type fastq --nproc $threads --sample_id ${sample} --bowtie2out "${metaphlan_dir}${sample}.bowtie2.bz2" "$f/stag-mwc/output_dir/host_removal/${sample}_1.fq.gz","$f/stag-mwc/output_dir/host_removal/${sample}_2.fq.gz" -o "${metaphlan_dir}${sample}.metaphlan.txt"   
    # looks like metaphlan 3 has broken downloads as well
    #
    # Convert MPA 3 output to something like MPA2 v2.7.7 output 
    # so it can be used with HUMAnN2, avoids StaG issue #138.
    # TODO: Remove this once HUMANn2 v2.9 is out.
    sed '/#/d' "${metaphlan_dir}${sample}.metaphlan.txt" | cut -f1,3 >> "${humann2_dir}mpa2_table-v2.7.7.txt"
    cat "$f/stag-mwc/output_dir/host_removal/${sample}_1.fq.gz" "$f/stag-mwc/output_dir/host_removal/${sample}_2.fq.gz" > "${humann2_dir}concat_input_reads.fq.gz"
    # humann2
    # ripping out humann2 to make the tests faster
    conda activate humann2
    humann2 --input "${humann2_dir}concat_input_reads.fq.gz" --output $humann2_dir --nucleotide-database "databases/func_databases/humann2/chocophlan" --protein-database "databases/func_databases/humann2/uniref" --output-basename $sample --threads $threads --taxonomic-profile "${humann2_dir}mpa2_table-v2.7.7.txt" 
    # normalize_humann2_tables
    humann2_renorm_table --input "${humann2_dir}${sample}_genefamilies.tsv" --output "${humann2_dir}${sample}_genefamilies_relab.tsv" --units relab --mode community 
    humann2_renorm_table --input "${humann2_dir}${sample}_pathabundance.tsv" --output "${humann2_dir}${sample}_pathabundance_relab.tsv" --units relab --mode community
    # join_humann2_tables
    humann2_join_tables --input $humann2_dir --output "${humann2_dir}all_samples.humann2_genefamilies.tsv" --file_name genefamilies_relab
    humann2_join_tables --input $humann2_dir --output "${humann2_dir}all_samples.humann2_pathabundance.tsv" --file_name pathcoverage
    humann2_join_tables --input $humann2_dir --output "${humann2_dir}all_samples.humann2_pathcoverage.tsv" --file_name pathabundance_relab
    # cleanup after finishing   
    rm "$f/stag-mwc/output_dir/humann2/concat_input_reads.fq.gz"
    rm -rf "$f/stag-mwc/output_dir/humann2/*_humann2_temp/" # the 1 isn't supposed to be static - it corresponds with the sample num
fi

mkdir -pv "$output_path/$name"

# Remove what is not needed for further analysis and takes up a lot of space
rm -rfv "$f/stag-mwc/output_dir/fastp/"
rm -rfv "$f/stag-mwc/output_dir/host_removal/"
rm -rfv "$f/stag-mwc/output_dir/amrplusplus/AlignToAMR/"
rm -fv "$f/stag-mwc/output_dir/kraken2/$name.kraken"

# Save the output folder and free up the space taken
datestamp=$(date -d "today" +"%Y%m%d%H%M")
mv "$f/stag-mwc/output_dir" "$output_path/$name/$datestamp"
cp "$f/stag-mwc/config.yaml" "$output_path/$name/$datestamp" # Config file might be useful for downstream analysis
for v in ENABLE_QC_READS ENABLE_HOST_REMOVAL ENABLE_KRAKEN2 ENABLE_GROOT ENABLE_AMRPLUSPLUS BRACKEN_TRESH ; do
    echo "$v=${!v}" >> "$output_path/$name/$datestamp/sub_params"
done
chmod g+w -R "$output_path/$name/$datestamp"
rm -rf "$f" # clean up after myself

echo "End time" $(date)
