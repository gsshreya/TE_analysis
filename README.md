# TE_analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21615045.svg)](https://doi.org/10.5281/zenodo.21615045)

Scripts and notebooks for our transposable element (TE) analysis, as described in the associated manuscript.

## MELT_Insertion_Calling/
MELT-based TE insertion calling pipeline, in order:
- `MELT_tosubmit.sh` : submits individual-sample calling jobs (uses `MELT_template.sh` per sample).
- `MELT_template.sh` : per-sample MELT Preprocess + IndivAnalysis template.
- `MoveIndivCallOutputOnly.sh` : moves individual-call output for a given sample list ahead of joint calling.
- `MELT_group_calling.sh` : MELT GroupAnalysis (joint calling) across all samples.
- `MELT_genotyping_inscreen.sh` : per-sample genotyping, run interactively in a `screen` session.
- `MELT_vcf_inscreen_alu.sh` : MELT MakeVCF for a given TE type (duplicated/edited per TE type).

## Benchmarking/
- `Benchmarking_Analysis.ipynb` : MELT vs. long-read (xTea+Sniffles) benchmarking: precision/recall/F1 by filter.

## Basic_Analysis/
- `InsertionandLengthAnalysis.ipynb` : Insertion counts, length distribution, AF distribution and PCA analysis.
- `Overlap_Analysis.ipynb` : Comparison and overlap with 8 other TE databases, as described in Table 1.
- `combined_pca.sh` : PCA on combined TE genotype data.

## Annotation/
- `te_annovar.sh` : ANNOVAR annotation of TE variants.
- `roadmap_rep.sh` : Roadmap Epigenomics annotation of TE insertions.
- `AnnotationsEnrichmentSummary.ipynb` : Summary of annotation and enrichment results.

## Enrichment/
- `enrichment_global.sh` : Fisher's exact enrichment of insertions in Roadmap chromatin-state annotations (promoters/enhancers/active/repressed, etc.) across tissues, ALU and LINE1.
- `ins_annovar_enr.sh` : Fisher's exact enrichment of insertions in ANNOVAR functional gene categories (exonic, intronic, intergenic, UTR3/5, splicing, ncRNA subtypes, upstream/downstream) vs. genome background, ALU and LINE1.

## Replication_Timing_Analysis/
- `replication_enrichment.sh` : Fisher's exact enrichment of common, non-singleton insertions in ERD/LRD replication-timing domains vs. genome background, pooled across all samples, ALU and LINE1.
- `pop_specific_repenr.sh` : same ERD/LRD enrichment analysis, run separately per population on non-singleton ALU insertions (one Fisher test per population group instead of pooled).
- `heterogeneity_test.R` : meta-analysis (`metafor`, fixed-effect Cochran's Q) comparing the TB group's Alu insertions ERD/LRD odds ratios against a pooled estimate across the other 4 groups, to test whether TB's replication-timing enrichment differs from the rest.

## GWAS_and_eQTL_Analysis/
- `ins_snps_ld.sh` : computes linkage disequilibrium between TE insertions and SNPs (pruned, merged VCF, PLINK LD).
- `GWAS_Associations_Final.ipynb` : Finding MEIs in LD with GWAS SNPs
- `eqtlscomparison_july2026.R` : Overlpping MEIs from this study with those found TE eQTLs found in GTEx
- `eQTL_Analysis_Final.ipynb` : Finding MEIs in LD with eQTLs and determining functional variants i.e. TE eQTLs, LD with eQTLs, LD with GWAS SNPs with differing AFs in 82 populations

## Exonic_TE/
- `exonic_TE_analysis.ipynb` : Annotates the functional region that each TE lies in and then filters out TE insertions falling in exonic regions. Filters exonic TE for further analysis.
- `CABYR_founder_analysis.py` : For the identified 6 carriers of insertion in gene CABYR, caculates haplotype level sharing in +/- 500 kb around the insertion breakpoint.


## Figures/
- 'fig_style.py': Style guides for figures generated for the manuscript.
- 'fig1.py': Used to generate figures 1A, 1B, 1C and 1D
- 'fig4_pca.py': Used to generate figures 4B and S9B

## Deletions/
- `melt_deletion.sh` : MELT indivudal deletion calling.
- `deletion_merge.sh` : MELT-Deletion-Merge script.
- `DeletionAnalysisMain.ipynb` : Deletion callset analysis and summary (deletion distrbution, benchmarking, AF distrbution)
