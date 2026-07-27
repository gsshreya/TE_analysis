#!/bin/bash
#SBATCH --job-name=MELT_jc
#SBATCH --partition=cbr_q_t
#SBATCH --exclusive
#SBATCH -c 72
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --output=/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/JointCall_7478/Group.%j.out
#SBATCH --error=/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/LOGS/Terminal_log/JointCall_7478/Group.%j.err


module load jdk1.8.0
module load bowtie2
module load samtools-1.14

#java -Xmx170G -jar /gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar GroupAnalysis -discoverydir /gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions -w /gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions -t /gpfs/data/user/shweta_lab/data/TE/TE_discovery/misc_files/transposon_ref_files.txt -h /gpfs/gidata/genomeindia/eric/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa -n /gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed
java -Xmx170G -jar /gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar GroupAnalysis \
-discoverydir /gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/ \
-w /gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/ \
-t /gpfs/data/user/shweta_lab/data/TE/TE_discovery/misc_files/transposon_ref_files.txt \
-h /gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-n /gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed