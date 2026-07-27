#!/bin/bash
#SBATCH --job-name=melt_deletion_genotype
#SBATCH -p cbr_q_huge
#SBATCH --array=0-2
#SBATCH -c 72
#SBATCH --mem=185G
#SBATCH -t 5-00:00:00
#SBATCH --exclusive
#SBATCH -o /gpfs/data/user/shreyags/TE_work/scripts/deletion/logs/melt_deletion_genotype_%A_%a.out
#SBATCH -e /gpfs/data/user/shreyags/TE_work/scripts/deletion/logs/melt_deletion_genotype_%A_%a.err

module load bcftools-1.21
module load samtools-1.21
module load parallel

source ~/.bashrc
conda activate rnaseq_r

LOGS="/gpfs/data/user/shreyags/TE_work/scripts/deletion/logs"
mkdir -p "$LOGS"

JAVA="$HOME/.conda/envs/rnaseq_r/bin/java"
MELT_JAR="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/MELT.jar"
CRAM_DIR="/gpfs/data/user/shweta_lab/data/bams/GI/markduplicates/crams"
REF="/gpfs/data/user/shweta_lab/data/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BED="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/tools/MELTv2.2.2/add_bed_files/Hg38/LINE1.AluY.deletion.bed"
OUTDIR="/gpfs/data/user/shreyags/TE_work/melt_deletion_f"
CHUNK_DIR="/gpfs/data/user/shreyags/TE_work/scripts/deletion"

mkdir -p "$OUTDIR"

THREADS=40
MEM_PER_JOB=4

SAMTOOLS=$(which samtools)
[ -z "$SAMTOOLS" ] && { echo "ERROR: samtools not found"; exit 1; }

FAIL_LOG="${OUTDIR}/failed_samples_${SLURM_ARRAY_TASK_ID}.txt"
touch "$FAIL_LOG"

echo "Using samtools: $SAMTOOLS"
echo "Array task: $SLURM_ARRAY_TASK_ID"

DONE_SAMPLES_FILE="${CHUNK_DIR}/done_samples_for_me.txt"

mapfile -t TASK_SAMPLES < "${CHUNK_DIR}/chunk_${SLURM_ARRAY_TASK_ID}.txt"
echo "This task processing ${#TASK_SAMPLES[@]} samples from chunk_${SLURM_ARRAY_TASK_ID}.txt"

TASK_CRAMS=()

for SAMPLE in "${TASK_SAMPLES[@]}"; do

    if grep -qx "$SAMPLE" "$DONE_SAMPLES_FILE"; then
        echo "Skipping (already done): $SAMPLE"
        continue
    fi

    CRAM="${CRAM_DIR}/${SAMPLE}.cram"
    if [ ! -f "$CRAM" ]; then
        echo "WARNING: CRAM not found for $SAMPLE — skipping"
        continue
    fi

    TASK_CRAMS+=("$CRAM")
done

run_deletion_genotype() {

    CRAM=$1
    MELT_JAR=$2
    BED=$3
    REF=$4
    MEM_PER_JOB=$5
    JAVA=$6
    SAMTOOLS=$7
    OUTDIR=$8

    SAMPLE=$(basename "$CRAM" .cram)
    ABS_CRAM=$(readlink -f "$CRAM")
    SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
    DEL_FILE="${SAMPLE_DIR}/${SAMPLE}.del.tsv"

    echo "Processing: $SAMPLE"

    cleanup_on_fail() {
        echo "Cleaning up partial output for $SAMPLE"
        rm -f "$DEL_FILE"
    }

    trap 'cleanup_on_fail' SIGINT SIGTERM EXIT

    if [ ! -f "$ABS_CRAM" ]; then
        echo "ERROR: CRAM missing $ABS_CRAM"
        echo -e "${SAMPLE}\tCRAM missing" >> "$FAIL_LOG"
        trap - EXIT
        return
    fi

    if [ ! -f "${ABS_CRAM}.crai" ]; then
        echo "ERROR: CRAM index missing ${ABS_CRAM}.crai"
        echo -e "${SAMPLE}\tCRAM index missing" >> "$FAIL_LOG"
        trap - EXIT
        return
    fi

    mkdir -p "$SAMPLE_DIR"

    LOCAL_CRAM="${SAMPLE_DIR}/${SAMPLE}.cram"
    ln -sf "$ABS_CRAM" "$LOCAL_CRAM"
    ln -sf "${ABS_CRAM}.crai" "${LOCAL_CRAM}.crai"

    $JAVA -Xmx${MEM_PER_JOB}g -jar "$MELT_JAR" Deletion-Genotype \
        -bamfile "$LOCAL_CRAM" \
        -bed "$BED" \
        -h "$REF" \
        -samtools "$SAMTOOLS" \
        -w "$SAMPLE_DIR" \
        2>> "${SAMPLE_DIR}/${SAMPLE}_deletion_genotype.log"

    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: Deletion-Genotype failed for $SAMPLE"
        echo -e "${SAMPLE}\tMELT failed (exit code $EXIT_CODE)" >> "$FAIL_LOG"
        cleanup_on_fail
    else
        echo "Done: $SAMPLE"
        trap - EXIT
    fi
}

export -f run_deletion_genotype
export JAVA MELT_JAR BED REF SAMTOOLS OUTDIR MEM_PER_JOB FAIL_LOG

parallel --line-buffer \
    --joblog "${LOGS}/parallel_joblog_${SLURM_ARRAY_TASK_ID}.log" \
    -j $THREADS run_deletion_genotype {} $MELT_JAR $BED $REF $MEM_PER_JOB $JAVA $SAMTOOLS $OUTDIR \
    ::: "${TASK_CRAMS[@]}"

echo "DELETION GENOTYPE COMPLETE for task $SLURM_ARRAY_TASK_ID"