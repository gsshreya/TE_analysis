#!/bin/bash

#Logging output and error

# Define log file paths
log_dir="/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/VCF"
mkdir -p "$log_dir"
output_log="${log_dir}/MELT_vcf_alu.out"
error_log="${log_dir}/MELT_vcf_alu.err"

# Redirect stdout and stderr to log files
exec > >(tee -a "$output_log")
exec 2> >(tee -a "$error_log" >&2)

module load jdk1.8.0
module load bowtie2
module load samtools-1.14
module load parallel

TE_name="ALU"
MELT_jar="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
ref_gen="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
work_dir="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/"
transposon_file_list="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/scripts/temp_lists_for_vcf_parallel/ALU_sublist.txt"
hg38_bed="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed"
#local ref_TE_bed="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/LINE1.AluY.deletion.bed"
# Perform processing here
java -jar -Xmx2G $MELT_jar MakeVCF \
-genotypingdir ${work_dir} \
-h ${ref_gen} \
-p ${work_dir} \
-t ${transposon_file_list} \
-w ${work_dir} 

echo "Processing $TE_name from $transposon_file_list"