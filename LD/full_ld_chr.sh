#!/bin/bash
#SBATCH --job-name=full_ld_chr
#SBATCH -p cbr_q_large
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 20:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

module load plink/2
module load plink/1.9.0


VCF="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/Tag_SNPs/Tag_SNPs_7478call/perchr/merged_autosomes.vcf.gz"

OUTDIR="$HOME/TE_work/replication_results/ld/per_chr" #results in /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Results/LD/per_chr
mkdir -p "$OUTDIR"
cd "$OUTDIR"

THREADS=${SLURM_CPUS_PER_TASK:-16}


if [[ -f merged.bed && -f merged.bim && -f merged.fam ]]; then
    echo "PLINK files already exist — skipping VCF conversion"
else
    echo "PLINK files not found — converting VCF → PLINK"

    plink2 \
      --vcf "$VCF" \
      --set-all-var-ids @:#:\$r:\$a \
      --make-bed \
      --allow-extra-chr \
      --threads $THREADS \
      --out merged
fi

#LD per chromosome

for CHR in {1..22}
do
    echo "=================================="
    echo "Chromosome $CHR"

    plink \
      --bfile merged \
      --chr $CHR \
      --r2 \
      --ld-window 100000000 \
      --ld-window-kb 100000000 \
      --ld-window-r2 0.7 \
      --threads $THREADS \
      --out ld_chr${CHR}
done

echo "Done!"