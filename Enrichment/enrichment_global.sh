#!/bin/bash
#SBATCH --job-name=TE_roadmap_tissue_fisher
#SBATCH -p cbr_q_small
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load bedtools-2.30
module load bcftools-1.21

source ~/.bashrc
conda activate rnaseq_r

mkdir -p logs

ALU_VCF="$HOME/TE_work/filtered_vcfs/ALU.filtered.final.vcf.gz"
LINE1_VCF="$HOME/TE_work/filtered_vcfs/LINE1.filtered.final.vcf.gz"

ROADMAP_DIR="/gpfs/data/user/shreyags/TE_work/data/roadmap_anno"
METADATA="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/EID_roadmap_metadata.tab"
GENOME="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/genome.bed"
OUTDIR="/gpfs/data/user/shreyags/TE_work/replication_results/new_enrichment/insertions/all"

mkdir -p "$OUTDIR"

genome_bp=$(awk '{s+=($3-$2)} END{print s+0}' "$GENOME")

declare -A TISSUES

while IFS=$'\t' read -r eid group color mnemonic std_name edacc_name anatomy type age sex solid ethnicity donor
do
    [[ "$eid" == "EID" ]] && continue
    TISSUES["$eid"]="$anatomy"
done < "$METADATA"

PROMOTER='^(1_TssA|2_TssAFlnk|10_TssBiv)$'
ENHANCER='^(6_EnhG|7_Enh|12_EnhBiv)$'
TRANSCRIBED='^(3_TxFlnk|4_Tx|5_TxWk)$'
REPRESSED='^(9_Het|13_ReprPC|14_ReprPCWk)$'
ACTIVE='^(1_TssA|2_TssAFlnk|3_TxFlnk|4_Tx|5_TxWk|6_EnhG|7_Enh|10_TssBiv|12_EnhBiv)$'

run_category() {

    NAME=$1
    REGEX=$2
    INS_BED=$3
    TE=$4
    TE_OUTDIR=$5

    COUNTS="${TE_OUTDIR}/${NAME}_counts.tsv"

    echo -e "tissue_id\ttissue_name\tfunctional_bp\tnonfunctional_bp\tinsertions_functional\tinsertions_nonfunctional\tnoninsertions_functional\tnoninsertions_nonfunctional\ttotal_insertions" > "$COUNTS"

    echo
    echo "${TE} : ${NAME}"

    for ANNO in ${ROADMAP_DIR}/E*_15_coreMarks_hg38lift_dense.bed.gz
    do

        tissue=$(basename "$ANNO" _15_coreMarks_hg38lift_dense.bed.gz)
        tissue_name="${TISSUES[$tissue]}"

        if [ -z "$tissue_name" ]; then
            tissue_name="UNKNOWN"
        fi

        BED="${TE_OUTDIR}/${tissue}.${NAME}.bed"

        echo "  ${tissue} (${tissue_name})"

        gzip -dc "$ANNO" |
        awk -v pat="$REGEX" '
            BEGIN{OFS="\t"}
            $1~/^chr/ && $4 ~ pat {
                print $1,$2,$3
            }
        ' |
        sort -k1,1 -k2,2n |
        bedtools merge -i - \
        > "$BED"

        functional_bp=$(awk '{s+=($3-$2)} END{print s+0}' "$BED")
        nonfunctional_bp=$((genome_bp - functional_bp))

        insertions_functional=$(bedtools intersect \
            -a "$INS_BED" \
            -b "$BED" \
            -u | wc -l)

        total_insertions=$(wc -l < "$INS_BED")
        insertions_nonfunctional=$((total_insertions - insertions_functional))
        noninsertions_functional=$((functional_bp - insertions_functional))
        noninsertions_nonfunctional=$((nonfunctional_bp - insertions_nonfunctional))

        echo -e "${tissue}\t${tissue_name}\t${functional_bp}\t${nonfunctional_bp}\t${insertions_functional}\t${insertions_nonfunctional}\t${noninsertions_functional}\t${noninsertions_nonfunctional}\t${total_insertions}" \
            >> "$COUNTS"

        rm -f "$BED"

    done

    Rscript - "$COUNTS" "${TE_OUTDIR}/${NAME}.tsv" << 'EOF'

args <- commandArgs(trailingOnly=TRUE)

infile  <- args[1]
outfile <- args[2]

x <- read.delim(infile, stringsAsFactors=FALSE)

for(col in c(
    "functional_bp",
    "nonfunctional_bp",
    "insertions_functional",
    "insertions_nonfunctional",
    "total_insertions"
)){
    x[[col]] <- as.numeric(x[[col]])
}

res_list <- vector("list", nrow(x))

SCALE_FACTOR <- 1000
THRESHOLD    <- 1e8

for(i in seq_len(nrow(x))) {

    row <- x[i,]

    A <- row$insertions_functional
    B <- row$total_insertions - row$insertions_functional

    C <- row$functional_bp
    D <- row$nonfunctional_bp

    # Only rescale C and D if they're large enough to risk exceeding
    # .Machine$integer.max (~2.1 billion)
    if (C > THRESHOLD || D > THRESHOLD) {
        C_scaled <- max(1, round(C / SCALE_FACTOR))
        D_scaled <- max(1, round(D / SCALE_FACTOR))
    } else {
        C_scaled <- C
        D_scaled <- D
    }

    mat <- matrix(
        c(A, B, C_scaled, D_scaled),
        nrow=2,
        byrow=TRUE
    )

    ft <- tryCatch(
        fisher.test(mat),
        error=function(e)
            list(
                estimate=NA,
                conf.int=c(NA,NA),
                p.value=NA
            )
    )

    res_list[[i]] <- data.frame(

        tissue_id   = row$tissue_id,
        tissue_name = row$tissue_name,

        functional_bp    = row$functional_bp,
        nonfunctional_bp = row$nonfunctional_bp,

        total_insertions = row$total_insertions,

        insertions_functional    = row$insertions_functional,
        insertions_nonfunctional = row$insertions_nonfunctional,

        A = A,
        B = B,
        C = C,
        D = D,

        C_scaled = C_scaled,
        D_scaled = D_scaled,

        prop_genome =
            row$functional_bp /
            (row$functional_bp + row$nonfunctional_bp),

        prop_insertions =
            row$insertions_functional /
            row$total_insertions,

        odds_ratio = unname(ft$estimate),
        ci_lower   = ft$conf.int[1],
        ci_upper   = ft$conf.int[2],
        pvalue     = ft$p.value,

        stringsAsFactors=FALSE
    )
}

results <- do.call(rbind, res_list)

results$FDR <- p.adjust(
    results$pvalue,
    method="BH"
)

results <- results[order(results$pvalue),]

write.table(
    results,
    outfile,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\nSaved:", outfile, "\n")

EOF

    rm -f "$COUNTS"

}

run_analysis() {

    TE=$1
    VCF=$2

    TE_OUTDIR="${OUTDIR}/${TE}"
    mkdir -p "$TE_OUTDIR"

    INS_BED="${TE_OUTDIR}/${TE}.insertions.bed"

    echo
    echo "Building insertion BED for ${TE}"

    bcftools query -f '%CHROM\t%POS\n' "$VCF" |
    awk 'BEGIN{OFS="\t"} $1~/^chr/ {print $1,$2-1,$2}' |
    sort -k1,1 -k2,2n \
    > "$INS_BED"

    run_category "promoters"   "$PROMOTER"    "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "enhancers"   "$ENHANCER"    "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "transcribed" "$TRANSCRIBED" "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "repressed"   "$REPRESSED"   "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "active"      "$ACTIVE"      "$INS_BED" "$TE" "$TE_OUTDIR"

    run_category "1_TssA"      '^1_TssA$'      "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "2_TssAFlnk"  '^2_TssAFlnk$'  "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "3_TxFlnk"    '^3_TxFlnk$'    "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "4_Tx"        '^4_Tx$'        "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "5_TxWk"      '^5_TxWk$'      "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "6_EnhG"      '^6_EnhG$'      "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "7_Enh"       '^7_Enh$'       "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "8_ZNF_Rpts"  '^8_ZNF/Rpts$'  "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "9_Het"       '^9_Het$'       "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "10_TssBiv"   '^10_TssBiv$'   "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "11_BivFlnk"  '^11_BivFlnk$'  "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "12_EnhBiv"   '^12_EnhBiv$'   "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "13_ReprPC"   '^13_ReprPC$'   "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "14_ReprPCWk" '^14_ReprPCWk$' "$INS_BED" "$TE" "$TE_OUTDIR"
    run_category "15_Quies"    '^15_Quies$'    "$INS_BED" "$TE" "$TE_OUTDIR"
}

run_analysis "ALU"   "$ALU_VCF"
run_analysis "LINE1" "$LINE1_VCF"

echo
echo "DONE"