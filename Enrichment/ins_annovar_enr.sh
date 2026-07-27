#!/bin/bash
#SBATCH --job-name=TE_annovar_fisher
#SBATCH -p cbr_q_small
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

mkdir -p logs

source ~/.bashrc
conda activate rnaseq_r

ALU_INS_ANNO="$HOME/TE_work/replication_results/insertions_annotations/ALU_ensGene_annovar.hg38_multianno.txt"
LINE1_INS_ANNO="$HOME/TE_work/replication_results/insertions_annotations/LINE1_ensGene_annovar.hg38_multianno.txt"

GENOME_CATEGORY_BP="$HOME/TE_work/replication_results/reference_annotations/genome_background/genome_category_bp.txt"
GENOME="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/genome.bed"

OUTDIR="$HOME/TE_work/replication_results/new_enrichment/annovar_ins_all"
mkdir -p "$OUTDIR"

FUNC_CATS=(
    "exonic"
    "intronic"
    "intergenic"
    "UTR3"
    "UTR5"
    "splicing"
    "ncRNA_exonic"
    "ncRNA_intronic"
    "ncRNA_splicing"
    "upstream"
    "downstream"
    "upstream:downstream"
)

RSCRIPT="${OUTDIR}/run_fisher.R"

cat > "$RSCRIPT" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]
outfile <- args[2]

x <- read.delim(infile, stringsAsFactors = FALSE)

for (col in c("variants_in","variants_out","background_in","background_out")) {
    x[[col]] <- as.numeric(x[[col]])
}

res_list <- vector("list", nrow(x))

SCALE_FACTOR <- 1000
THRESHOLD    <- 1e8

for (i in seq_len(nrow(x))) {

    row <- x[i, ]

    A <- row$variants_in
    B <- row$variants_out
    C <- row$background_in
    D <- row$background_out

    # Only rescale C and D for insertion comparisons (vs_genome), where
    # background values are genome-scale bp counts that can exceed
    # .Machine$integer.max (~2.1 billion).
    if (row$comparison == "vs_genome" && (C > THRESHOLD || D > THRESHOLD)) {
        C_scaled <- max(1, round(C / SCALE_FACTOR))
        D_scaled <- max(1, round(D / SCALE_FACTOR))
    } else {
        C_scaled <- C
        D_scaled <- D
    }

    mat <- matrix(
        as.numeric(c(A, B, C_scaled, D_scaled)),
        nrow  = 2,
        byrow = TRUE
    )

    ft <- tryCatch(
        fisher.test(mat),
        error = function(e) list(
            estimate = NA,
            conf.int = c(NA, NA),
            p.value  = NA
        )
    )

    total_variants   <- A + B
    total_background <- C + D

    res_list[[i]] <- data.frame(
        te_type          = row$te_type,
        category         = row$category,
        comparison       = row$comparison,
        total_variants   = total_variants,
        total_background = total_background,
        variants_in      = A,
        background_in    = C,
        background_in_scaled  = C_scaled,
        background_out_scaled = D_scaled,
        prop_variants    = ifelse(total_variants   > 0, A / total_variants,   NA),
        prop_background  = ifelse(total_background > 0, C / total_background, NA),
        odds_ratio       = unname(ft$estimate),
        ci_lower         = ft$conf.int[1],
        ci_upper         = ft$conf.int[2],
        pvalue           = ft$p.value,
        stringsAsFactors = FALSE
    )
}

results     <- do.call(rbind, res_list)
results$FDR <- p.adjust(results$pvalue, method = "BH")

write.table(results, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:", outfile, "\n")
REOF

get_col_index() {
    local FILE=$1
    local COLNAME=$2
    head -1 "$FILE" | tr '\t' '\n' | grep -n "^${COLNAME}$" | cut -d: -f1
}

count_annovar_cats() {
    local MULTIANNO=$1
    local OUTFILE=$2
    local FUNC_COL=$3

    awk -v fc="$FUNC_COL" '
        NR == 1 { next }
        {
            n = split($fc, parts, ";")
            for (p = 1; p <= n; p++) {
                cat = parts[p]
                gsub(/^ +| +$/, "", cat)
                if (cat != "" && cat != ".") func_count[cat]++
            }
        }
        END { for (c in func_count) print c "\t" func_count[c] }
    ' "$MULTIANNO" | sort > "$OUTFILE"
}

append_category_row() {
    local CAT=$1
    local TE_TYPE=$2
    local COMPARISON=$3
    local VAR_COUNTS=$4
    local TOTAL_VARIANTS=$5
    local BG_COUNTS=$6
    local TOTAL_BG=$7
    local COUNTS_FILE=$8

    awk -v cat="$CAT" \
        -v te="$TE_TYPE" \
        -v comp="$COMPARISON" \
        -v vcounts="$VAR_COUNTS" \
        -v total_var="$TOTAL_VARIANTS" \
        -v bcounts="$BG_COUNTS" \
        -v total_bg="$TOTAL_BG" \
    'BEGIN {
        while ((getline line < vcounts) > 0) {
            split(line, a, "\t")
            vmap[a[1]] = a[2]+0
        }
        while ((getline line < bcounts) > 0) {
            split(line, a, "\t")
            bmap[a[1]] = a[2]+0
        }
        var_in  = (cat in vmap) ? vmap[cat] : 0
        var_out = total_var - var_in
        bg_in   = (cat in bmap) ? bmap[cat] : 0
        bg_out  = total_bg - bg_in
        print te "\t" cat "\t" comp "\t" var_in "\t" var_out "\t" bg_in "\t" bg_out
    }' >> "$COUNTS_FILE"
}

run_insertion() {
    local TE=$1
    local MULTIANNO=$2
    local COUNTS_FILE=$3

    echo "INSERTIONS ${TE}"

    local FUNC_COL
    FUNC_COL=$(get_col_index "$MULTIANNO" "Func.ensGene")
    if [[ -z "$FUNC_COL" ]]; then
        echo "ERROR: Func.ensGene not found in $MULTIANNO"; return 1
    fi

    local VAR_COUNTS
    VAR_COUNTS=$(mktemp)
    count_annovar_cats "$MULTIANNO" "$VAR_COUNTS" "$FUNC_COL"

    local TOTAL_VARIANTS
    TOTAL_VARIANTS=$(awk 'NR>1' "$MULTIANNO" | wc -l)

    local GENOME_BP
    GENOME_BP=$(awk '{s+=($3-$2)} END{print s+0}' "$GENOME")

    for CAT in "${FUNC_CATS[@]}"; do
        append_category_row \
            "$CAT" "${TE}_ins" "vs_genome" \
            "$VAR_COUNTS" "$TOTAL_VARIANTS" \
            "$GENOME_CATEGORY_BP" "$GENOME_BP" \
            "$COUNTS_FILE"
    done

    rm -f "$VAR_COUNTS"
}

COUNTS_FILE="${OUTDIR}/insertions_counts.tsv"
echo -e "te_type\tcategory\tcomparison\tvariants_in\tvariants_out\tbackground_in\tbackground_out" \
    > "$COUNTS_FILE"

run_insertion "ALU"   "$ALU_INS_ANNO"   "$COUNTS_FILE"
run_insertion "LINE1" "$LINE1_INS_ANNO" "$COUNTS_FILE"

echo
echo "Running Fisher tests..."

RESULT_FILE="${OUTDIR}/insertions_results.tsv"
Rscript "$RSCRIPT" "$COUNTS_FILE" "$RESULT_FILE"

for CAT in "${FUNC_CATS[@]}"; do
    SAFE_CAT="${CAT//\//_}"
    PER_CAT_FILE="${OUTDIR}/${SAFE_CAT}_results.tsv"
    awk -F'\t' -v cat="$CAT" '
        NR == 1 { print; next }
        $2 == cat { print }
    ' "$RESULT_FILE" > "$PER_CAT_FILE"
done

COMPILED_CSV="${OUTDIR}/compiled_summary.csv"

awk -F'\t' 'BEGIN{OFS=","}
    NR == 1 {
        for (i=1; i<=NF; i++) col[$i] = i
        print "te_type","category","odds_ratio","ci_lower","ci_upper","pvalue"
        next
    }
    {
        print $col["te_type"], $col["category"], $col["odds_ratio"], \
              $col["ci_lower"], $col["ci_upper"], $col["pvalue"]
    }
' "$RESULT_FILE" > "$COMPILED_CSV"

echo "DONE"
