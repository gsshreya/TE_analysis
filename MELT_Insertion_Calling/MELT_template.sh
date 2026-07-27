#!/bin/bash
#SBATCH --job-name=MELT
#SBATCH --partition=cbr_q_t
#SBATCH --exclusive
#SBATCH -c 72
#SBATCH --nodes=1
#SBATCH --mem=187G
#SBATCH --output=/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/57_rerun_7480/MELT.%j.out
#SBATCH --error=/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/57_rerun_7480/MELT.%j.err

error="/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/MELT_log/57_rerun_7480/"
MELT_jar="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
ref_gen="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
work_dir="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/"

transposon_file_list="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/misc_files/transposon_ref_files.txt"
hg38_bed="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed"

module load jdk1.8.0
module load bowtie2
module load samtools-1.14

sample_list="{{SAMPLES}}"

for sample in $sample_list; do
  echo "Processing sample: ${sample}"
  
  # Run MELT commands in the background with logging for each sample
  ( 
    java -jar -Xmx8G $MELT_jar Preprocess \
        -bamfile "${work_dir}${sample}" \
        -h "${ref_gen}" && \
    java -jar -Xmx8G $MELT_jar IndivAnalysis \
        -w "${work_dir}" \
        -bamfile "${work_dir}${sample}" \
        -t "${transposon_file_list}" \
        -h "${ref_gen}" && \
    rm -f "${work_dir}${sample}.disc" "${work_dir}${sample}.disc.bai" "${work_dir}${sample}.fq"
  ) > "${error}${sample}.out" 2> "${error}${sample}.err" &
done

# Wait for all background processes to complete
wait