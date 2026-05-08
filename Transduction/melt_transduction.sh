#!/bin/bash
#SBATCH --job-name=melt_transduction_find
#SBATCH -p cbr_q_large
#SBATCH --array=0-1
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 72:00:00
#SBATCH -o logs/melt_transduction_%A_%a.out
#SBATCH -e logs/melt_transduction_%A_%a.err

module load bcftools-1.21
module load samtools-1.9
module load parallel
source ~/.bashrc
conda activate rnaseq_r

mkdir -p logs

JAVA="$HOME/.conda/envs/rnaseq_r/bin/java"
MELT_JAR="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
CRAM_DIR="/gpfs/data/user/shweta_lab/data/bams/GI/markduplicates/crams"
REF="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SOURCE_BED="$HOME/TE_work/source_beds/hg38.filtered.FL.canonical.bed.source"
JOINT_VCF="$HOME/TE_work/filtered_vcfs/LINE1.filtered.vcf.gz"
OUTDIR="$HOME/TE_work/melt_transduction_f"
INDEX_DIR="$OUTDIR/cram_indices"   # where we store indices we create ourselves

mkdir -p "$OUTDIR" "$INDEX_DIR"

THREADS=14
MEM_PER_JOB=8

SAMTOOLS=$(which samtools)
if [ -z "$SAMTOOLS" ]; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

echo "Using samtools: $SAMTOOLS"
echo "Array task: $SLURM_ARRAY_TASK_ID"

# Build full valid CRAM list
VCF_SAMPLES=$(bcftools query -l $JOINT_VCF)
VALID_CRAMS=()

for CRAM in $CRAM_DIR/*.cram; do
    SAMPLE=$(basename $CRAM .cram)

    if ! echo "$VCF_SAMPLES" | grep -qx "$SAMPLE"; then
        continue
    fi

    ABS_CRAM=$(readlink -f "$CRAM")

    # Check if index exists next to cram
    if [ ! -f "${ABS_CRAM}.crai" ] && [ ! -f "${ABS_CRAM%.cram}.crai" ]; then
        # No index found so create it in INDEX_DIR
        CRAI_PATH="${INDEX_DIR}/${SAMPLE}.cram.crai"
        if [ ! -f "$CRAI_PATH" ]; then
            echo "No CRAI found for $SAMPLE — indexing into $INDEX_DIR ..."
            $SAMTOOLS index -@ 4 "$ABS_CRAM" "$CRAI_PATH"
            if [ $? -ne 0 ]; then
                echo "ERROR: Failed to index $SAMPLE — skipping"
                continue
            fi
        else
            echo "Using pre-built index in INDEX_DIR for $SAMPLE"
        fi
    fi

    VALID_CRAMS+=("$CRAM")
done

TOTAL=${#VALID_CRAMS[@]}
echo "Total valid samples: $TOTAL"

# Split into 2 halves based on array task ID
HALF=$(( (TOTAL + 1) / 2 ))
START=$(( SLURM_ARRAY_TASK_ID * HALF ))
END=$(( START + HALF ))
if [ $END -gt $TOTAL ]; then END=$TOTAL; fi

TASK_CRAMS=("${VALID_CRAMS[@]:$START:$((END - START))}")
echo "This task processing samples $START to $END (${#TASK_CRAMS[@]} samples)"

run_transduction_find() {
    CRAM=$1
    MELT_JAR=$2
    SOURCE_BED=$3
    REF=$4
    MEM_PER_JOB=$5
    JAVA=$6
    SAMTOOLS=$7
    OUTDIR=$8
    INDEX_DIR=$9

    SAMPLE=$(basename $CRAM .cram)
    ABS_CRAM=$(readlink -f "$CRAM")
    SAMPLE_DIR="${OUTDIR}/${SAMPLE}"

    # Skip if already processed
    if [ -f "${SAMPLE_DIR}/${SAMPLE}.trans" ]; then
        echo "${SAMPLE} .trans file found, moving to next"
        return
    fi

    echo "Processing: $SAMPLE"

    if [ ! -f "$ABS_CRAM" ]; then
        echo "ERROR: CRAM missing $ABS_CRAM"
        return
    fi

    mkdir -p "$SAMPLE_DIR"

    # Symlink CRAM
    LOCAL_CRAM="${SAMPLE_DIR}/${SAMPLE}.cram"
    ln -sf "$ABS_CRAM" "$LOCAL_CRAM"

    # Symlink index prefer one next to original CRAM, fall back to INDEX_DIR
    if [ -f "${ABS_CRAM}.crai" ]; then
        ln -sf "${ABS_CRAM}.crai" "${LOCAL_CRAM}.crai"
    elif [ -f "${ABS_CRAM%.cram}.crai" ]; then
        ln -sf "${ABS_CRAM%.cram}.crai" "${LOCAL_CRAM}.crai"
    elif [ -f "${INDEX_DIR}/${SAMPLE}.cram.crai" ]; then
        ln -sf "${INDEX_DIR}/${SAMPLE}.cram.crai" "${LOCAL_CRAM}.crai"
    else
        echo "ERROR: No CRAI found for $SAMPLE even after indexing step — skipping"
        return
    fi

    $JAVA -Xmx${MEM_PER_JOB}g -jar $MELT_JAR TransductionFind \
        -bamfile "$LOCAL_CRAM" \
        -source "$SOURCE_BED" \
        -h "$REF" \
        -samtools "$SAMTOOLS" \
        2> "${SAMPLE_DIR}/${SAMPLE}_find.log"

    if [ $? -ne 0 ]; then
        echo "ERROR: TransductionFind failed for $SAMPLE"
    else
        echo "Done: $SAMPLE"
    fi
}

export -f run_transduction_find

parallel -j $THREADS run_transduction_find {} $MELT_JAR $SOURCE_BED $REF $MEM_PER_JOB $JAVA $SAMTOOLS $OUTDIR $INDEX_DIR ::: "${TASK_CRAMS[@]}"

echo "=== STEP 1 COMPLETE for task $SLURM_ARRAY_TASK_ID ==="