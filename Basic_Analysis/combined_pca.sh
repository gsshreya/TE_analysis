#!/bin/bash
#SBATCH --job-name=TE_PCA
#SBATCH -p cbr_q_small
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

source ~/.bashrc
conda activate rnaseq_r

module load plink/2
module load htslib-1.18

#VCF=/gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Filtered_VCFs/ALU_LINE1.filtered.tagged.vcf.gz
#VCF=/gpfs/data/user/shreyags/TE_work/filtered_vcfs/common_nosingletons_new_final/ALU_LINE1.filtered.final.no_singletons.vcf.gz
VCF=/gpfs/data/user/shreyags/TE_work/filtered_vcfs/with_af_new_final/ALU_LINE1.filtered.final.tagged.vcf.gz
OUTDIR=$HOME/TE_work/analysis/INS_PCA_all
POPMAP=/gpfs/data/user/shreyags/TE_work/filtered_vcfs/sample_population_map.csv

PLINK_PREFIX="${OUTDIR}/ALU_LINE1.merged"
PRUNE_PREFIX="${OUTDIR}/ALU_LINE1.merged.pruned"
PCA_PREFIX="${OUTDIR}/ALU_LINE1.merged.pca"

# LD pruning parameters (--indep-pairwise <window_kb> <step> <r2_threshold>)
LD_WINDOW_KB=50
LD_STEP=5
LD_R2=0.2

# Number of PCs to compute
N_PCS=10

mkdir -p "${OUTDIR}" logs

module load bcftools-1.18

echo "[$(date '+%F %T')] Starting PLINK PCA job"
echo "  Input VCF   : ${VCF}"
echo "  OUTDIR      : ${OUTDIR}"
echo "  Pop map     : ${POPMAP}"
echo "  LD pruning  : window=${LD_WINDOW_KB}kb step=${LD_STEP} r2=${LD_R2}"
echo "  No MAF filter applied (using all variants)"

if [[ ! -f "$VCF" ]]; then
    echo "[FATAL] VCF not found: $VCF"
    exit 1
fi

if [[ ! -f "$POPMAP" ]]; then
    echo "[FATAL] Population map not found: $POPMAP"
    exit 1
fi

FILTERED_VCF="${OUTDIR}/ALU_LINE1.merged.no_chrposrefalt_dups.vcf.gz"

echo "[$(date '+%F %T')] Removing duplicate CHR:POS:REF:ALT records..."

bcftools view -h "${VCF}" > "${OUTDIR}/vcf_header.tmp"

bcftools view -H "${VCF}" | \
awk '
BEGIN { OFS="\t" }
{
    key = $1 ":" $2 ":" $4 ":" $5
    if (!(key in seen)) {
        seen[key] = 1
        print
    }
}
' > "${OUTDIR}/vcf_body.tmp"

cat "${OUTDIR}/vcf_header.tmp" "${OUTDIR}/vcf_body.tmp" | \
bgzip -c > "${FILTERED_VCF}"

tabix -f -p vcf "${FILTERED_VCF}"

rm "${OUTDIR}/vcf_header.tmp" "${OUTDIR}/vcf_body.tmp"

echo "[$(date '+%F %T')] Importing filtered VCF into PLINK2..."

plink2 \
    --vcf "${FILTERED_VCF}" \
    --allow-extra-chr \
    --set-all-var-ids '@:#:$r:$a' \
    --vcf-half-call m \
    --make-pgen \
    --out "${PLINK_PREFIX}"

if [[ ! -f "${PLINK_PREFIX}.pgen" ]]; then
    echo "[FATAL] PLINK2 import failed -- no .pgen produced. Check ${PLINK_PREFIX}.log"
    exit 1
fi

N_VARIANTS_IN=$(wc -l < "${PLINK_PREFIX}.pvar")
echo "[$(date '+%F %T')] Imported variants (incl. header line): ${N_VARIANTS_IN}"

#LD-based pruning

echo "[$(date '+%F %T')] Running LD pruning (--indep-pairwise ${LD_WINDOW_KB} ${LD_STEP} ${LD_R2})..."

plink2 \
    --pfile "${PLINK_PREFIX}" \
    --allow-extra-chr \
    --indep-pairwise "${LD_WINDOW_KB}" "${LD_STEP}" "${LD_R2}" \
    --out "${PRUNE_PREFIX}"

if [[ ! -f "${PRUNE_PREFIX}.prune.in" ]]; then
    echo "[FATAL] LD pruning failed -- no .prune.in produced. Check ${PRUNE_PREFIX}.log"
    exit 1
fi

N_KEPT=$(wc -l < "${PRUNE_PREFIX}.prune.in")
N_REMOVED=$(wc -l < "${PRUNE_PREFIX}.prune.out" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] Variants kept after pruning : ${N_KEPT}"
echo "[$(date '+%F %T')] Variants removed by pruning  : ${N_REMOVED}"

if [[ "$N_KEPT" -eq 0 ]]; then
    echo "[FATAL] No variants survived LD pruning."
    exit 1
fi

#PCA on the pruned variant set

echo "[$(date '+%F %T')] Running PCA on pruned variant set (${N_PCS} PCs)..."

plink2 \
    --pfile "${PLINK_PREFIX}" \
    --allow-extra-chr \
    --extract "${PRUNE_PREFIX}.prune.in" \
    --pca "${N_PCS}" \
    --out "${PCA_PREFIX}"

if [[ ! -f "${PCA_PREFIX}.eigenvec" ]]; then
    echo "[FATAL] PCA failed -- no .eigenvec produced. Check ${PCA_PREFIX}.log"
    exit 1
fi

echo "[$(date '+%F %T')] PCA complete."
echo "  ${PCA_PREFIX}.eigenvec  -- per-sample PC scores"
echo "  ${PCA_PREFIX}.eigenval  -- eigenvalues (variance per PC)"


cat > postprocess_pca.R << 'EOF'

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: postprocess_pca.R <pca_prefix> <output_dir> <population_map.csv>")
}

PCA_PREFIX  <- args[1]
OUTDIR      <- args[2]
POPMAP_PATH <- args[3]

eigenvec_path <- paste0(PCA_PREFIX, ".eigenvec")
eigenval_path <- paste0(PCA_PREFIX, ".eigenval")

cat(sprintf("[%s] Reading eigenvalues: %s\n", ts(), eigenval_path))
eigenval <- scan(eigenval_path, what = numeric())

var_explained <- eigenval / sum(eigenval)
n_pcs <- length(var_explained)

var_table <- data.frame(
  PC = paste0("PC", seq_len(n_pcs)),
  eigenvalue = eigenval,
  variance_explained = var_explained,
  cumulative_variance = cumsum(var_explained)
)
write.csv(var_table,
          file.path(OUTDIR, "pca_variance_explained.csv"),
          row.names = FALSE)
cat(sprintf("[%s] Wrote variance-explained table (%d PCs)\n", ts(), n_pcs))


cat(sprintf("[%s] Reading sample PC scores: %s\n", ts(), eigenvec_path))
scores <- read.table(eigenvec_path, header = TRUE, comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE)

id_col <- grep("IID$", colnames(scores), value = TRUE)[1]
if (is.na(id_col)) stop("Could not find an IID column in eigenvec file.")
colnames(scores)[colnames(scores) == id_col] <- "sample_id"

pc_cols <- grep("^PC[0-9]+$", colnames(scores), value = TRUE)
scores <- scores[, c("sample_id", pc_cols)]


cat(sprintf("[%s] Reading population map: %s\n", ts(), POPMAP_PATH))
popmap <- read.csv(POPMAP_PATH, header = TRUE, stringsAsFactors = FALSE,
                    colClasses = "character")

if (ncol(popmap) < 4) {
  stop(sprintf("Population map has only %d columns -- expected at least 4 (sample, ?, population, treatment-flag).",
               ncol(popmap)))
}

map_sample <- trimws(popmap[[1]])
map_col2   <- trimws(popmap[[2]])
map_pop    <- trimws(popmap[[3]])
map_flag   <- trimws(popmap[[4]])

build_label <- function(pop, flag) {

  if (is.na(flag) || is.na(pop)) {
    if (is.na(pop)) return(NA_character_)
    return(pop)
  }

  flag_norm <- toupper(flag)
  if (flag_norm %in% c("N/A", "NA", "")) {
    return(pop)
  } else if (flag_norm == "YES") {
    return(paste0(pop, "_T"))
  } else if (flag_norm == "NO") {
    return(paste0(pop, "_NT"))
  } else {
    return(paste0(pop, "_UNK[", flag, "]"))
  }
}

map_label <- mapply(build_label, map_pop, map_flag)
popmap_lookup_group1 <- setNames(map_label, map_sample)
popmap_lookup_group2 <- setNames(map_col2, map_sample)

strip_sample_name <- function(s) sub("_.*$", "", s)

stripped_ids <- strip_sample_name(scores$sample_id)
scores$group1_pop_treatment <- popmap_lookup_group1[stripped_ids]
scores$group2_col2          <- popmap_lookup_group2[stripped_ids]

report_unmatched <- function(values, label) {
  n_unmatched <- sum(is.na(values))
  if (n_unmatched > 0) {
    cat(sprintf("[%s] [WARN] %d / %d samples did not match the population map for %s grouping (stripped IDs shown):\n",
                ts(), n_unmatched, length(values), label))
    unmatched_ids <- unique(stripped_ids[is.na(values)])
    cat(sprintf("         %s\n", paste(head(unmatched_ids, 20), collapse = ", ")))
    if (length(unmatched_ids) > 20) {
      cat(sprintf("         ... and %d more\n", length(unmatched_ids) - 20))
    }
    values[is.na(values)] <- "UNMATCHED"
  }
  values
}

scores$group1_pop_treatment <- report_unmatched(scores$group1_pop_treatment, "col3_T/NT")
scores$group2_col2          <- report_unmatched(scores$group2_col2, "col2")

write.csv(scores,
          file.path(OUTDIR, "pca_sample_scores.csv"),
          row.names = FALSE)
cat(sprintf("[%s] Wrote sample PC scores with both groupings: group1_pop_treatment (col3_T/NT) and group2_col2 (%d PCs, %d samples)\n",
            ts(), length(pc_cols), nrow(scores)))

cat(sprintf("[%s] Done. Outputs written to: %s\n", ts(), OUTDIR))
cat(sprintf("  - pca_sample_scores.csv       (per-sample PC1-PC%d + both groupings -- INPUT for plot_pca.R)\n", length(pc_cols)))
cat(sprintf("  - pca_variance_explained.csv  (eigenvalue + %% variance per PC -- INPUT for plot_pca.R)\n"))
EOF

Rscript postprocess_pca.R "${PCA_PREFIX}" "${OUTDIR}" "${POPMAP}"

echo "[$(date '+%F %T')] TE_PCA job complete. Results in ${OUTDIR}"