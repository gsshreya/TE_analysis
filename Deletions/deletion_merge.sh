#!/bin/bash
#SBATCH --job-name=melt_deletion_merge
#SBATCH -p cbr_q_t
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/melt_deletion_merge_%j.out
#SBATCH -e logs/melt_deletion_merge_%j.err

module load samtools-1.21

source ~/.bashrc
conda activate rnaseq_r

mkdir -p logs

JAVA="$HOME/.conda/envs/rnaseq_r/bin/java"
MELT_JAR="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"

REF="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BED="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/LINE1.AluY.deletion.bed"

OUTDIR="$HOME/TE_work/melt_deletion_f"
MERGE_DIR="$HOME/TE_work/melt_deletion_merge_new"

mkdir -p "$MERGE_DIR"

DEL_LIST="${MERGE_DIR}/deletion_files.txt"

find "$OUTDIR" -name "*.del.tsv" | sort > "$DEL_LIST"

TOTAL=$(wc -l < "$DEL_LIST")

echo "Found $TOTAL deletion TSV files"

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No .del.tsv files found"
    exit 1
fi

echo "Running MELT DeletionMerge..."

$JAVA -Xmx24g -jar "$MELT_JAR" Deletion-Merge \
    -bed "$BED" \
    -h "$REF" \
    -mergelist "$DEL_LIST" \
    -o "$MERGE_DIR"

EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: MELT DeletionMerge failed"
    exit $EXIT_CODE
fi

echo "DELETION MERGE COMPLETE"