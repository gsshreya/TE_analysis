#!/bin/bash

# Script to process BAM files using MELT without SLURM directives

# Set variables
sample_list="/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/scripts/temp_lists_for_genotyping_parallel/ToBeReprocessed_batch003.txt"
num_samples_parallel=50
total_samples_inthisJOB=271

# Define log file paths
log_dir="/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/Genotyping"
mkdir -p "$log_dir"
output_log="${log_dir}/MELT_genotype_batch003_Reprocess.out"
error_log="${log_dir}/MELT_genotype_batch003_Reprocess.err"

# Redirect stdout and stderr to log files
exec > >(tee -a "$output_log")
exec 2> >(tee -a "$error_log" >&2)

echo "Processing $num_samples_parallel samples parallel from $sample_list with $total_samples_inthisJOB"

# Load necessary modules
module load jdk1.8.0
module load bowtie2
module load samtools-1.14
module load parallel

# Define the function to process BAM files
process_bam() {
    local bam_file="$1"
    local MELT_jar="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
    local ref_gen="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    local work_dir="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/"
    local transposon_file_list="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/misc_files/transposon_ref_files.txt"
    local hg38_bed="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed"
    #local ref_TE_bed="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/LINE1.AluY.deletion.bed"
    sleep 10s
    if [ ! -e "${work_dir}${bam_file%.bam.bai}.ALU.tsv" ]; then
        # Perform processing
        java -jar -Xmx2G $MELT_jar Genotype -bamfile ${work_dir}${bam_file%.bai} -w ${work_dir} -p ${work_dir} -h ${ref_gen} -t ${transposon_file_list}
        echo "Processing ${bam_file%.bam.bai}"
        sleep 10s
    fi
}

# Export the function for parallel execution
export -f process_bam

# Run the function in parallel for the provided sample list
parallel -j $num_samples_parallel process_bam :::: "$sample_list"