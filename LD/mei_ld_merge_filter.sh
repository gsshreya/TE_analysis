#!/bin/bash
#SBATCH --job-name=mei_ld_filter_merge
#SBATCH -p cbr_q_small
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load bcftools-1.21

LINE1_VCF="$HOME/TE_work/filtered_vcfs/LINE1.filtered.vcf.gz" #these files are also in /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Filtered_VCFS
ALU_VCF="$HOME/TE_work/filtered_vcfs/ALU.filtered.vcf.gz" #these files are also in /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Filtered_VCFS
SVA_VCF="/gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/scripts/final_vcfs_7478/SVA.final_comp.vcf.gz"

LD_DIR="$HOME/TE_work/replication_results/ld/per_chr" #equivalent to /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Results/LD/per_chr
OUTDIR="$HOME/TE_work/replication_results/ld" #results in /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Results/LD/

mkdir -p $OUTDIR


bcftools query -f '%CHROM:%POS:%REF:%ALT\n' $LINE1_VCF \
| sed 's/^chr//' \
> $OUTDIR/LINE1_meis.txt


bcftools query -f '%CHROM:%POS:%REF:%ALT\n' $ALU_VCF \
| sed 's/^chr//' \
> $OUTDIR/ALU_meis.txt


bcftools query -f '%CHROM:%POS:%REF:%ALT\t%FILTER\t%INFO/ASSESS\n' $SVA_VCF \
| awk '$2=="PASS" && $3==5 && $1 !~ /_/' \
| cut -f1 \
| sed 's/^chr//' \
> $OUTDIR/SVA_meis.txt


cat $OUTDIR/*_meis.txt | sort -u > $OUTDIR/valid_meis.txt


OUTLD="$OUTDIR/ld_filtered_merged.ld"
> $OUTLD


for LD in $LD_DIR/*.ld
do
    echo "Processing $LD"

    awk 'NR==FNR{valid[$1]; next}
  {
      meiA = ($3 ~ /INS:ME:/)
      meiB = ($6 ~ /INS:ME:/)

      A = ($3 in valid)
      B = ($6 in valid)

      if(!meiA && !meiB) {print; next}
      if(meiA && !meiB && A) {print; next}
      if(!meiA && meiB && B) {print; next}
      if(meiA && meiB && A && B) {print}
  }
  ' $OUTDIR/valid_meis.txt $LD >> $OUTLD

done


echo "Merged filtered LD pairs:"
wc -l $OUTLD