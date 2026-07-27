#!/bin/bash
#SBATCH --job-name=reptime_fisher
#SBATCH -p cbr_q_small
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load bedtools-2.30
module load bcftools-1.21

source ~/.bashrc
conda activate rnaseq_r

mkdir -p logs

#ALU_INS_VCF=~/TE_work/filtered_vcfs/ALU.filtered.final.vcf.gz
ALU_INS_VCF="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/common_nosingletons_new_final/ALU.filtered.final.common_maf0.01.vcf.gz"

#LINE1_INS_VCF=~/TE_work/filtered_vcfs/LINE1.filtered.final.vcf.gz
LINE1_INS_VCF="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/common_nosingletons_new_final/LINE1.filtered.final.common_maf0.01.vcf.gz"

DOMAINS=~/TE_work/replication_results/enrichment/replication/Bg02es_segments_hg38.bed

GENOME="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/genome.bed"

OUT_DIR=~/TE_work/replication_results/new_enrichment/replication_common
mkdir -p "$OUT_DIR"

echo "Extracting ERD and LRD..."

grep "ERD" "$DOMAINS" | sort -k1,1 -k2,2n > "${OUT_DIR}/ERD_hg38.bed"
grep "LRD" "$DOMAINS" | sort -k1,1 -k2,2n > "${OUT_DIR}/LRD_hg38.bed"

ERD_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${OUT_DIR}/ERD_hg38.bed")
LRD_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${OUT_DIR}/LRD_hg38.bed")
GENOME_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$GENOME")

echo "ERD_BP=$ERD_BP  LRD_BP=$LRD_BP  GENOME_BP=$GENOME_BP"

#Insertion BEDs

echo "Converting insertion VCFs to BED..."

for TYPE in ALU LINE1; do
    if [[ "$TYPE" == "ALU" ]]; then
        VCF="$ALU_INS_VCF"
    else
        VCF="$LINE1_INS_VCF"
    fi

    bcftools query -f '%CHROM\t%POS\n' "$VCF" |
    awk 'BEGIN{OFS="\t"} $1~/^chr/ {print $1, $2-1, $2}' \
    > "${OUT_DIR}/${TYPE}_insertion.bed"

    echo "  ${TYPE} insertions: $(wc -l < ${OUT_DIR}/${TYPE}_insertion.bed)"
done

#Intersect variants with ERD/LRD

echo "Intersecting variants with ERD/LRD..."

for TYPE in ALU LINE1; do
    BED="${OUT_DIR}/${TYPE}_insertion.bed"
    bedtools intersect -a "$BED" -b "${OUT_DIR}/ERD_hg38.bed" -u > "${OUT_DIR}/${TYPE}_insertion_in_ERD.bed"
    bedtools intersect -a "$BED" -b "${OUT_DIR}/LRD_hg38.bed" -u > "${OUT_DIR}/${TYPE}_insertion_in_LRD.bed"
    echo "  ${TYPE} insertion in ERD: $(wc -l < ${OUT_DIR}/${TYPE}_insertion_in_ERD.bed)"
    echo "  ${TYPE} insertion in LRD: $(wc -l < ${OUT_DIR}/${TYPE}_insertion_in_LRD.bed)"
done

# Insertions: background = genome bp
#   2x2: variants vs genome bp, in ERD/LRD vs not

COUNTS="${OUT_DIR}/fisher_counts.tsv"
echo -e "TE\tClass\tVar_ERD\tVar_LRD\tVar_Total\tBg_ERD\tBg_LRD\tBg_Total\tBg_Type" > "$COUNTS"

for TYPE in ALU LINE1; do

    VAR_ERD=$(wc -l < "${OUT_DIR}/${TYPE}_insertion_in_ERD.bed")
    VAR_LRD=$(wc -l < "${OUT_DIR}/${TYPE}_insertion_in_LRD.bed")
    VAR_TOTAL=$(wc -l < "${OUT_DIR}/${TYPE}_insertion.bed")

    echo -e "${TYPE}\tinsertion\t${VAR_ERD}\t${VAR_LRD}\t${VAR_TOTAL}\t${ERD_BP}\t${LRD_BP}\t${GENOME_BP}\tbp" >> "$COUNTS"

done

Rscript - "$COUNTS" "${OUT_DIR}/fisher_results.tsv" <<'EOF'

args <- commandArgs(trailingOnly = TRUE)
dat  <- read.delim(args[1], stringsAsFactors = FALSE)
out_file <- args[2]

for (col in c("Var_ERD","Var_LRD","Var_Total","Bg_ERD","Bg_LRD","Bg_Total")) {
    dat[[col]] <- as.numeric(dat[[col]])
}

run_test <- function(A, B, C, D) {
    mat <- matrix(as.numeric(c(A, B, C, D)), nrow = 2, byrow = TRUE)
    ft  <- tryCatch(
        fisher.test(mat),
        error = function(e) list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
    )
    c(or = unname(ft$estimate), ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2], p = ft$p.value)
}

res <- vector("list", nrow(dat))

for (i in seq_len(nrow(dat))) {

    row <- dat[i, ]

    # ERD fisher
    # Rows: variants | background
    # Cols: in_ERD   | not_in_ERD
    A_erd <- row$Var_ERD
    B_erd <- row$Var_Total - row$Var_ERD
    C_erd <- row$Bg_ERD
    D_erd <- row$Bg_Total - row$Bg_ERD
    erd   <- run_test(A_erd, B_erd, C_erd, D_erd)

    # LRD fisher
    A_lrd <- row$Var_LRD
    B_lrd <- row$Var_Total - row$Var_LRD
    C_lrd <- row$Bg_LRD
    D_lrd <- row$Bg_Total - row$Bg_LRD
    lrd   <- run_test(A_lrd, B_lrd, C_lrd, D_lrd)

    res[[i]] <- data.frame(
        TE              = row$TE,
        Class           = row$Class,
        Background      = row$Bg_Type,
        Var_Total       = row$Var_Total,
        Bg_Total        = row$Bg_Total,
        ERD_Var         = A_erd,
        ERD_Bg          = C_erd,
        ERD_Prop_Var    = A_erd / row$Var_Total,
        ERD_Prop_Bg     = C_erd / row$Bg_Total,
        ERD_OR          = erd["or"],
        ERD_CI_lo       = erd["ci_lo"],
        ERD_CI_hi       = erd["ci_hi"],
        ERD_P           = erd["p"],
        LRD_Var         = A_lrd,
        LRD_Bg          = C_lrd,
        LRD_Prop_Var    = A_lrd / row$Var_Total,
        LRD_Prop_Bg     = C_lrd / row$Bg_Total,
        LRD_OR          = lrd["or"],
        LRD_CI_lo       = lrd["ci_lo"],
        LRD_CI_hi       = lrd["ci_hi"],
        LRD_P           = lrd["p"],
        stringsAsFactors = FALSE
    )
}

results          <- do.call(rbind, res)
results$ERD_FDR  <- p.adjust(results$ERD_P, method = "BH")
results$LRD_FDR  <- p.adjust(results$LRD_P, method = "BH")

write.table(results, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:", out_file, "\n")
print(results)

EOF

echo "Done. Results in ${OUT_DIR}/fisher_results.tsv"