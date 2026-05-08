#!/bin/bash
#SBATCH --job-name=mei_pca
#SBATCH -p cbr_q_small
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

source ~/.bashrc
conda activate rnaseq_r

module load bcftools-1.21

INPUTDIR="/gpfs/data/user/shreyags/TE_work/filtered_vcfs"
LINE1_VCF="$INPUTDIR/LINE1.filtered.vcf.gz"
ALU_VCF="$INPUTDIR/ALU.filtered.vcf.gz"
SAMPLE_MAP="$INPUTDIR/sample_population_map.csv"
OUTDIR="/gpfs/data/user/shreyags/TE_work/replication_results/pca"

mkdir -p "$OUTDIR" logs

bcftools query -H -f '%CHROM:%POS[\t%GT]\n' "$LINE1_VCF" > "$OUTDIR/LINE1.gt.tsv"
bcftools query -H -f '%CHROM:%POS[\t%GT]\n' "$ALU_VCF"   > "$OUTDIR/ALU.gt.tsv"

export OUTDIR
export SAMPLE_MAP

Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

run_mei_pca <- function(gt_file, sample_map) {

  gt <- read_tsv(gt_file, show_col_types = FALSE)

  # strip bcftools -H prefix "# [1]" from first col, "[N]" from rest
  cn <- colnames(gt)
  cn <- sub("^#\\s*", "", cn)
  cn <- sub("^\\[[0-9]+\\]", "", cn)
  cn <- sub(":GT$", "", cn)
  cn[-1] <- sub("_.*$", "", cn[-1])
  colnames(gt) <- cn

  mat <- as.matrix(gt[, -1])
  rownames(mat) <- gt[[1]]

  # handle both phased and unphased genotypes
  mat[mat %in% c("0/0", "0|0")] <- 0
  mat[mat %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  mat[mat %in% c("1/1", "1|1")] <- 2
  mat[mat %in% c("./.", ".|.")] <- NA

  storage.mode(mat) <- "numeric"
  mat <- t(mat)   # samples x variants

  # impute missing with column mean
  col_means <- colMeans(mat, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) {
    na_rows <- is.na(mat[, j])
    if (any(na_rows)) mat[na_rows, j] <- col_means[j]
  }

  # drop zero-variance sites
  mat <- mat[, apply(mat, 2, sd) > 0, drop = FALSE]

  # --- MAF FILTER (comment out this block to skip) ---
  #maf_threshold <- 0.01
  #af  <- colMeans(mat, na.rm = TRUE) / 2
  #maf <- pmin(af, 1 - af)
  #mat <- mat[, maf > maf_threshold, drop = FALSE]
  #cat(sprintf("  Variants after MAF>%.2f filter: %d\n", maf_threshold, ncol(mat)))
  # --- END MAF FILTER ---

  pca     <- prcomp(mat, center = TRUE, scale. = TRUE)
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

  # read metadata, index by name to satisfy dplyr/.data pronoun rules
  raw <- read_csv(sample_map, show_col_types = FALSE)
  nms <- colnames(raw)

  sm <- raw |>
    transmute(
      sample       = sub("_.*$", "", as.character(.data[[nms[1]]])),
      population   = as.character(.data[[nms[2]]]),
      superpop_raw = as.character(.data[[nms[3]]]),
      twin_flag    = as.character(.data[[nms[4]]])
    ) |>
    filter(!is.na(sample), superpop_raw != "") |>
    mutate(
      super_population = paste0(
        superpop_raw, "_",
        ifelse(twin_flag == "Yes", "T", "NT")
      )
    ) |>
    distinct(sample, .keep_all = TRUE)

  df <- as.data.frame(pca$x[, 1:4]) |>
    rownames_to_column("sample") |>
    left_join(sm, by = "sample") |>
    filter(!is.na(population))

  df$population       <- factor(df$population,       levels = sort(unique(df$population)))
  df$super_population <- factor(df$super_population, levels = sort(unique(df$super_population)))

  list(pca = pca, df = df, var = var_exp)
}

plot_pca_grid <- function(obj, colour, title, outfile) {

  var    <- obj$var
  groups <- levels(obj$df[[colour]])

  pc_pairs <- list(
    c("PC1","PC2"), c("PC2","PC3"),
    c("PC1","PC3"), c("PC3","PC4"),
    c("PC2","PC4"), c("PC1","PC4")
  )

  pal <- setNames(
    scales::hue_pal(l = 65, c = 100)(length(groups)),
    groups
  )

  plots <- lapply(pc_pairs, function(pcs) {
    pcx <- as.integer(sub("PC", "", pcs[1]))
    pcy <- as.integer(sub("PC", "", pcs[2]))

    ggplot(obj$df,
           aes(.data[[pcs[1]]], .data[[pcs[2]]],
               color = .data[[colour]])) +
      geom_point(size = 1.6, alpha = 0.9) +
      scale_color_manual(values = pal, drop = FALSE) +
      theme_bw() +
      labs(
        x     = paste0(pcs[1], " (", round(var[pcx] * 100, 2), "%)"),
        y     = paste0(pcs[2], " (", round(var[pcy] * 100, 2), "%)"),
        color = colour
      )
  })

  # number of legend columns — scales with number of groups
  legend_ncol <- ceiling(length(groups) / 10)

  final_plot <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 6),
      legend.key.size  = unit(0.3, "cm"),
      legend.title     = element_text(size = 7)
    )

  final_plot <- final_plot +
    plot_annotation(title = title) &
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 2)))

  ggsave(outfile, final_plot, width = 16, height = 14)
}

OUTDIR     <- Sys.getenv("OUTDIR")
SAMPLE_MAP <- Sys.getenv("SAMPLE_MAP")

pca_LINE1 <- run_mei_pca(file.path(OUTDIR, "LINE1.gt.tsv"), SAMPLE_MAP)
pca_ALU   <- run_mei_pca(file.path(OUTDIR, "ALU.gt.tsv"),   SAMPLE_MAP)

plot_pca_grid(pca_LINE1, "population",       "LINE1 Population PCA", file.path(OUTDIR, "LINE1_population.pdf"))
plot_pca_grid(pca_ALU,   "population",       "ALU Population PCA",   file.path(OUTDIR, "ALU_population.pdf"))
plot_pca_grid(pca_LINE1, "super_population", "LINE1 Superpop PCA",   file.path(OUTDIR, "LINE1_superpop.pdf"))
plot_pca_grid(pca_ALU,   "super_population", "ALU Superpop PCA",     file.path(OUTDIR, "ALU_superpop.pdf"))

cat("DONE\n")
EOF