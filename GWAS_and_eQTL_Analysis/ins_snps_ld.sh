#!/bin/bash
#SBATCH --job-name=ins_snps_ld
#SBATCH -p cbr_q_small
#SBATCH --nodelist node9
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load bcftools-1.18
module load plink/2
module load jdk-20.0.2
module load bedtools-2.30
module load plink/1.9.0
module load htslib-1.18

set -euo pipefail

Alu_ins_vcf=~/TE_work/filtered_vcfs/common_nosingletons_new_final/ALU.filtered.final.common_maf0.01.vcf.gz
L1_ins_vcf=~/TE_work/filtered_vcfs/common_nosingletons_new_final/LINE1.filtered.final.common_maf0.01.vcf.gz

work_dir=/tmp/${USER}_insertions_snps_LD/

# Final destination for results on gpfs
RESULTS_DIR=~/TE_work/analysis/insertions_snps_LD_results/

SNP_vcf=/gpfs/data/user/shweta_lab/data/GI/9772callset/concatenatedautosomes_9772.vcf.gz

mkdir -p ${work_dir}temp
mkdir -p ${RESULTS_DIR}

bcftools query -l ${Alu_ins_vcf} | cut -f 1 -d '_' | sort | uniq > ${work_dir}temp/TE_samples.txt
bcftools query -l ${SNP_vcf} | sort | uniq > ${work_dir}temp/SNP_samples.txt

grep -F -f ${work_dir}temp/TE_samples.txt ${work_dir}temp/SNP_samples.txt | sort | uniq > ${work_dir}samples_list.txt && \
rm -f ${work_dir}temp/TE_samples.txt ${work_dir}temp/SNP_samples.txt && samples_list=${work_dir}samples_list.txt

# Subset insertion VCFs to shared samples and set variant IDs
bcftools query -l $Alu_ins_vcf | cut -f 1 -d '_' | bcftools reheader -s /dev/stdin $Alu_ins_vcf | \
    bcftools annotate -h <(echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">') | \
    bcftools view -S ${samples_list} | \
    bcftools annotate --set-id +'%CHROM\:%POS\:%FIRST_ALT' -Oz \
    -o ${work_dir}temp/Alu_ins_common.vcf.gz
tabix -p vcf ${work_dir}temp/Alu_ins_common.vcf.gz

bcftools query -l $L1_ins_vcf | cut -f 1 -d '_' | bcftools reheader -s /dev/stdin $L1_ins_vcf | \
    bcftools annotate -h <(echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">') | \
    bcftools view -S ${samples_list} | \
    bcftools annotate --set-id +'%CHROM\:%POS\:%FIRST_ALT' -Oz \
    -o ${work_dir}temp/L1_ins_common.vcf.gz
tabix -p vcf ${work_dir}temp/L1_ins_common.vcf.gz

# Create ±100kb windows around insertions for SNP extraction
zcat ${work_dir}temp/Alu_ins_common.vcf.gz | grep -v '#' | awk 'BEGIN{OFS="\t"}{
    start = ($2 - 100001 < 0) ? 0 : $2 - 100001
    print $1, start, $2+100000
}' > ${work_dir}temp/alu_ins_pos.bed

zcat ${work_dir}temp/L1_ins_common.vcf.gz | grep -v '#' | awk 'BEGIN{OFS="\t"}{
    start = ($2 - 100001 < 0) ? 0 : $2 - 100001
    print $1, start, $2+100000
}' > ${work_dir}temp/l1_ins_pos.bed

# Merge overlapping windows
cat ${work_dir}temp/alu_ins_pos.bed ${work_dir}temp/l1_ins_pos.bed | sort -k1,1 -k2,2n | bedtools merge -i - \
    > ${work_dir}regions.bed && \
    rm -f ${work_dir}temp/alu_ins_pos.bed ${work_dir}temp/l1_ins_pos.bed

echo "regions.bed covers $(awk '{s+=$3-$2} END{print s}' ${work_dir}regions.bed) bp"

# Extract SNPs in those regions for shared samples, filter to MAF > 0.01
bcftools view -R ${work_dir}regions.bed -S ${work_dir}samples_list.txt ${SNP_vcf} | \
    bcftools +fill-tags -- -t MAF | \
    bcftools view -i 'MAF>0.01' | \
    bcftools annotate --set-id +'%CHROM\:%POS\:%FIRST_ALT' | \
    bgzip -@ 8 > ${work_dir}SNP.vcf.gz && \
    tabix -p vcf ${work_dir}SNP.vcf.gz

# Prune SNPs for LD
plink2 --vcf ${work_dir}SNP.vcf.gz --mind 0.1 --geno 0.1 --indep-pairwise 50 5 0.5 --out ${work_dir}temp/snp --threads 8
plink2 --vcf ${work_dir}SNP.vcf.gz --extract ${work_dir}temp/snp.prune.in --make-bed --out ${work_dir}temp/snp_pruned --threads 8
plink2 --bfile ${work_dir}temp/snp_pruned --export vcf bgz --out ${work_dir}snp_pruned --threads 8

rm -f ${work_dir}SNP.vcf.gz ${work_dir}SNP.vcf.gz.tbi
rm -f ${work_dir}temp/snp_pruned.bed ${work_dir}temp/snp_pruned.bim ${work_dir}temp/snp_pruned.fam

plink2 --vcf ${work_dir}temp/Alu_ins_common.vcf.gz --geno 0.1 --make-bed --out ${work_dir}temp/temp --threads 8
plink2 --bfile ${work_dir}temp/temp --export vcf bgz --out ${work_dir}Alu_ins_pruned

plink2 --vcf ${work_dir}temp/L1_ins_common.vcf.gz --geno 0.1 --make-bed --out ${work_dir}temp/temp --threads 8
plink2 --bfile ${work_dir}temp/temp --export vcf bgz --out ${work_dir}L1_ins_pruned

rm -f ${work_dir}temp/temp.bed ${work_dir}temp/temp.bim ${work_dir}temp/temp.fam
rm -f ${work_dir}temp/Alu_ins_common.vcf.gz* ${work_dir}temp/L1_ins_common.vcf.gz*

# Reheader pruned SNP VCF to match insertion VCF header
bcftools view -h ${work_dir}Alu_ins_pruned.vcf.gz | \
    bcftools reheader -h /dev/stdin ${work_dir}snp_pruned.vcf.gz | \
    bcftools sort -Oz -o ${work_dir}snp_pruned1.vcf.gz

rm -f ${work_dir}snp_pruned.vcf.gz

java -jar /gpfs/data/user/shweta_lab/bin/picard.jar MergeVcfs \
    -I ${work_dir}Alu_ins_pruned.vcf.gz \
    -I ${work_dir}L1_ins_pruned.vcf.gz \
    -I ${work_dir}snp_pruned1.vcf.gz \
    -O ${work_dir}merged.vcf.gz \
    --CREATE_INDEX false

tabix -p vcf ${work_dir}merged.vcf.gz

MERGED_SIZE_BYTES=$(stat -c%s ${work_dir}merged.vcf.gz)
MERGED_SIZE_GB=$((MERGED_SIZE_BYTES / 1024 / 1024 / 1024))
echo "merged.vcf.gz (pruned MEI insertions + pruned SNPs) size: ${MERGED_SIZE_GB} GB (${MERGED_SIZE_BYTES} bytes)"

# Compute LD, report pairs with r2 >= 0.7
plink2 --vcf ${work_dir}merged.vcf.gz --make-bed --out ${work_dir}temp/for_ld --threads 8
plink --bfile ${work_dir}temp/for_ld \
    --r2 \
    --ld-window 100000000 \
    --ld-window-kb 100000000 \
    --ld-window-r2 0.7 \
    --out ${work_dir}ld_results

cp ${work_dir}ld_results.* ${RESULTS_DIR}

echo "Done. Only ld_results.* copied to ${RESULTS_DIR}"
echo "All intermediate/merged files remain in ${work_dir} only until node cleanup"