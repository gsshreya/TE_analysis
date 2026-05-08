#!/bin/bash
#SBATCH --job-name=transduction_merge
#SBATCH -p cbr_q_large
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 72:00:00
#SBATCH -o logs/transduction_merge_%j.out
#SBATCH -e logs/transduction_merge_%j.err

module load bcftools-1.21
module load samtools-1.21
source ~/.bashrc
conda activate rnaseq_r

mkdir -p logs

JAVA="$HOME/.conda/envs/rnaseq_r/bin/java"
MELT_JAR="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
CRAM_DIR="/gpfs/data/user/shweta_lab/data/bams/GI/markduplicates/crams"
REF="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SOURCE_BED="$HOME/TE_work/source_beds/hg38.filtered.FL.canonical.bed.source"
JOINT_VCF="$HOME/TE_work/filtered_vcfs/LINE1.filtered.7338.vcf.gz"
OUTDIR="$HOME/TE_work/melt_transduction_f"
BAMLIST="$OUTDIR/bamlist_for_merge.txt"
COVERAGE_DIR="/gpfs/data/user/ruchitha/scripts/CoverageCalculation/temperoray_files_coverage_mosdepth"

SAMTOOLS=$(which samtools)
if [ -z "$SAMTOOLS" ]; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

echo "Using samtools: $SAMTOOLS"

# Build bamlist
VCF_SAMPLES=$(bcftools query -l $JOINT_VCF)
> "$BAMLIST"

SKIPPED_NO_TRANS=0
SKIPPED_NO_CRAM=0
SKIPPED_NO_COV=0

for CRAM in $CRAM_DIR/*.cram; do
    SAMPLE=$(basename $CRAM .cram)

    # Must be in VCF
    if ! echo "$VCF_SAMPLES" | grep -qx "$SAMPLE"; then
        continue
    fi

    ABS_CRAM=$(readlink -f "$CRAM")
    SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
    LOCAL_CRAM="${SAMPLE_DIR}/${SAMPLE}.cram"

    # Check symlinked CRAM exists in sample dir
    if [ ! -f "$LOCAL_CRAM" ]; then
        echo "WARNING: No CRAM found at $LOCAL_CRAM — skipping $SAMPLE"
        (( SKIPPED_NO_CRAM++ ))
        continue
    fi

    # Check .trans file exists in sample dir
    if [ ! -f "${SAMPLE_DIR}/${SAMPLE}.trans" ]; then
        echo "WARNING: No .trans file for $SAMPLE — skipping"
        (( SKIPPED_NO_TRANS++ ))
        continue
    fi

    # Strip _md suffix for coverage file lookup
    SAMPLE_BASE="${SAMPLE/_md/}"
    COV_FILE="${COVERAGE_DIR}/${SAMPLE_BASE}_coverage.txt"
    if [ ! -f "$COV_FILE" ]; then
        echo "WARNING: No coverage file for $SAMPLE (looked for ${SAMPLE_BASE}_coverage.txt) — skipping"
        (( SKIPPED_NO_COV++ ))
        continue
    fi

    COV=$(grep "^${SAMPLE_BASE}" "$COV_FILE" | awk '{gsub("X","",$2); print int($2)}')
    if [ -z "$COV" ] || [ "$COV" -eq 0 ]; then
        echo "WARNING: Could not parse coverage for $SAMPLE — skipping"
        (( SKIPPED_NO_COV++ ))
        continue
    fi

    echo -e "${LOCAL_CRAM}\t${COV}" >> "$BAMLIST"
    echo "Added $SAMPLE (coverage: ${COV}x)"
done

TOTAL=$(wc -l < "$BAMLIST")
echo ""
echo "=== Bamlist summary ==="
echo "Samples added:         $TOTAL"
echo "Skipped (no CRAM):     $SKIPPED_NO_CRAM"
echo "Skipped (no .trans):   $SKIPPED_NO_TRANS"
echo "Skipped (no coverage): $SKIPPED_NO_COV"
echo "Bamlist written to:    $BAMLIST"
echo ""

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No valid samples in bamlist — exiting"
    exit 1
fi

# Decompress VCF 
VCF_INPUT="$JOINT_VCF"
if [[ "$JOINT_VCF" == *.gz ]]; then
    VCF_PLAIN="${OUTDIR}/$(basename ${JOINT_VCF%.gz})"
    if [ ! -f "$VCF_PLAIN" ]; then
        echo "Decompressing VCF for MELT input..."
        bcftools view "$JOINT_VCF" > "$VCF_PLAIN"
    fi
    VCF_INPUT="$VCF_PLAIN"
fi

# Run TransductionMerge 
echo "Running TransductionMerge on $TOTAL samples..."

$JAVA -Xmx12G -jar $MELT_JAR TransductionMerge \
    -bamlist "$BAMLIST" \
    -h "$REF" \
    -source "$SOURCE_BED" \
    -vcf "$VCF_INPUT"

if [ $? -ne 0 ]; then
    echo "ERROR: TransductionMerge failed"
    exit 1
fi

echo "=== TransductionMerge complete ==="
echo "Output: ${OUTDIR}/LINE1.filtered.trans.vcf"