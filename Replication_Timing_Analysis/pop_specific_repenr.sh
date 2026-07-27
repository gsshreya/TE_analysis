#!/bin/bash
#SBATCH --job-name=alu_ins_popbypop_reptime
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

VCF_DIR="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/pops/5groups_nosingletons"

DOMAINS="$HOME/TE_work/replication_results/enrichment/replication/Bg02es_segments_hg38.bed"

GENOME="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/genome.bed"

OUT_DIR="$HOME/TE_work/replication_results/new_enrichment/pop_replication_all_nosingletons"

mkdir -p "$OUT_DIR"

echo "Extracting ERD and LRD..."

grep "ERD" "$DOMAINS" | sort -k1,1 -k2,2n > "${OUT_DIR}/ERD_hg38.bed"
grep "LRD" "$DOMAINS" | sort -k1,1 -k2,2n > "${OUT_DIR}/LRD_hg38.bed"

ERD_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${OUT_DIR}/ERD_hg38.bed")
LRD_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${OUT_DIR}/LRD_hg38.bed")
GENOME_BP=$(awk '{sum += $3-$2} END {print sum+0}' "$GENOME")

echo "ERD bp: $ERD_BP"
echo "LRD bp: $LRD_BP"
echo "Genome bp: $GENOME_BP"

COUNTS="${OUT_DIR}/counts.tsv"

echo -e \
"Dataset\tVar_ERD\tVar_LRD\tVar_Total\tBg_ERD\tBg_LRD\tBg_Total" \
> "$COUNTS"

echo
echo "Processing VCFs..."

# Loop over each pop subdir and find its ALU_INS VCF
for POP_DIR in "${VCF_DIR}"/*/
do
    POP=$(basename "$POP_DIR")

    VCF="${POP_DIR}${POP}_ALU_INS.vcf.gz"

    [[ ! -f "$VCF" ]] && echo "  SKIPPING $POP -- VCF not found" && continue

    DATASET="$POP"

    echo "  $DATASET"

    BED="${OUT_DIR}/${DATASET}.bed"

    bcftools query \
        -f '%CHROM\t%POS\n' \
        "$VCF" |
    awk 'BEGIN{OFS="\t"}
        $1~/^chr/{
            print $1,$2-1,$2
        }' \
    > "$BED"

    bedtools intersect \
        -a "$BED" \
        -b "${OUT_DIR}/ERD_hg38.bed" \
        -u \
    > "${OUT_DIR}/${DATASET}_ERD.bed"

    bedtools intersect \
        -a "$BED" \
        -b "${OUT_DIR}/LRD_hg38.bed" \
        -u \
    > "${OUT_DIR}/${DATASET}_LRD.bed"

    VAR_ERD=$(wc -l < "${OUT_DIR}/${DATASET}_ERD.bed")
    VAR_LRD=$(wc -l < "${OUT_DIR}/${DATASET}_LRD.bed")
    VAR_TOTAL=$(wc -l < "$BED")

    echo -e \
"${DATASET}\t${VAR_ERD}\t${VAR_LRD}\t${VAR_TOTAL}\t${ERD_BP}\t${LRD_BP}\t${GENOME_BP}" \
    >> "$COUNTS"

done

echo
echo "Running Fisher tests..."

Rscript - "$COUNTS" "${OUT_DIR}/subset_replication_results.tsv" << 'EOF'

args <- commandArgs(trailingOnly=TRUE)

counts_file <- args[1]
out_file <- args[2]

dat <- read.delim(
    counts_file,
    stringsAsFactors = FALSE
)

for(col in c(
    "Var_ERD",
    "Var_LRD",
    "Var_Total",
    "Bg_ERD",
    "Bg_LRD",
    "Bg_Total"
)){
    dat[[col]] <- as.numeric(dat[[col]])
}

run_test <- function(A,B,C,D){

    ft <- fisher.test(
        matrix(
            c(A,B,C,D),
            nrow=2,
            byrow=TRUE
        )
    )

    c(
        A = A, B = B, C = C, D = D,
        OR    = unname(ft$estimate),
        CI_LO = ft$conf.int[1],
        CI_HI = ft$conf.int[2],
        P     = ft$p.value
    )
}

results_list <- vector(
    "list",
    nrow(dat)
)

for(i in seq_len(nrow(dat))){

    row <- dat[i,]

    erd <- run_test(
        row$Var_ERD,
        row$Var_Total - row$Var_ERD,
        row$Bg_ERD,
        row$Bg_Total - row$Bg_ERD
    )

    lrd <- run_test(
        row$Var_LRD,
        row$Var_Total - row$Var_LRD,
        row$Bg_LRD,
        row$Bg_Total - row$Bg_LRD
    )

    results_list[[i]] <- data.frame(

        Dataset = row$Dataset,

        Var_Total = row$Var_Total,

        #ERD 2x2 cell counts 
        ERD_Var_In     = erd["A"],   # insertions in ERD
        ERD_Var_Out    = erd["B"],   # insertions NOT in ERD
        ERD_Bg_In      = erd["C"],   # background bp in ERD
        ERD_Bg_Out     = erd["D"],   # background bp NOT in ERD

        ERD_Count = row$Var_ERD,
        ERD_Prop_Var = row$Var_ERD / row$Var_Total,
        ERD_Prop_Bg = row$Bg_ERD / row$Bg_Total,
        ERD_OR = erd["OR"],
        ERD_CI_lo = erd["CI_LO"],
        ERD_CI_hi = erd["CI_HI"],
        ERD_P = erd["P"],

        #LRD 2x2 cell counts 
        LRD_Var_In     = lrd["A"],   # insertions in LRD
        LRD_Var_Out    = lrd["B"],   # insertions NOT in LRD
        LRD_Bg_In      = lrd["C"],   # background bp in LRD
        LRD_Bg_Out     = lrd["D"],   # background bp NOT in LRD

        LRD_Count = row$Var_LRD,
        LRD_Prop_Var = row$Var_LRD / row$Var_Total,
        LRD_Prop_Bg = row$Bg_LRD / row$Bg_Total,
        LRD_OR = lrd["OR"],
        LRD_CI_lo = lrd["CI_LO"],
        LRD_CI_hi = lrd["CI_HI"],
        LRD_P = lrd["P"],

        stringsAsFactors = FALSE

    )

}

results <- do.call(
    rbind,
    results_list
)

results$ERD_FDR <- p.adjust(
    results$ERD_P,
    method="BH"
)

results$LRD_FDR <- p.adjust(
    results$LRD_P,
    method="BH"
)

results <- results[
    order(results$Dataset),
]

write.table(
    results,
    out_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\nSaved:", out_file, "\n")

print(results)

EOF

echo "DONE"

echo
echo "Results:"
echo "${OUT_DIR}/subset_replication_results.tsv"