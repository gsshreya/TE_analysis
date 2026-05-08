#!/bin/bash
#SBATCH --job-name=GI_Roadmap
#SBATCH --partition=cbr_q_large
#SBATCH --exclusive
#SBATCH -c 72
#SBATCH --nodes=1
#SBATCH --mem=187G
#SBATCH --output=logs/GI_Roadmap.%j.out
#SBATCH --error=logs/GI_Roadmap.%j.err

module load bedtools-2.30
module load parallel
module load bcftools-1.21

GI_VCF_DIR="$HOME/TE_work/filtered_vcfs"
ROADMAP_PATH="/gpfs/data/user/shreyags/TE_work/data/roadmap_anno"
OUTDIR="$HOME/TE_work/replication_results/annotations/roadmap/GI"

mkdir -p "$OUTDIR"

# MEI TYPES TO PROCESS

MEIS=("ALU" "LINE1")

# ROADMAP FUNCTION

annotate_roadmap () {

    local roadmapfile=$1
    local GI_BED=$2
    local MEI=$3

    name=$(basename "$roadmapfile")
    tissue=$(echo "$name" | sed 's/_15_coreMarks_hg38lift_dense.bed.gz//')

    outfile="${OUTDIR}/GI_${MEI}_${tissue}_Roadmap.tsv.gz"

    bedtools intersect \
        -a "$GI_BED" \
        -b <(zcat "$roadmapfile") \
        -wa -wb \
    | awk 'BEGIN{OFS="\t"}{
        print $1 ":" ($2+1), $7
    }' \
    | gzip -c > "$outfile"

    echo "Saved $outfile"
}

export -f annotate_roadmap
export OUTDIR


# MAIN LOOP


for MEI in "${MEIS[@]}"; do

  echo "=== Processing GI ${MEI} ==="

  TMPDIR="/gpfs/data/user/$USER/tmp/GI_Roadmap_${MEI}_${SLURM_JOB_ID}"
  mkdir -p "$TMPDIR"

  GI_VCF="${GI_VCF_DIR}/${MEI}.filtered.vcf.gz"
  GI_BED="${TMPDIR}/GI_${MEI}.bed"

  echo "Preparing BED for GI ${MEI}"

  bcftools view -H "$GI_VCF" \
  | awk 'BEGIN{OFS="\t"}{
      chr=$1;
      if(chr !~ /^chr/) chr="chr"chr;
      print chr, $2-1, $2
  }' > "$GI_BED"

  find "$ROADMAP_PATH" -name "*_15_coreMarks_hg38lift_dense.bed.gz" \
  | parallel -j 72 annotate_roadmap {} "$GI_BED" "$MEI"

  rm -rf "$TMPDIR"

  echo "Completed GI ${MEI}"

done

echo "Roadmap annotation complete for GI ALU and LINE1"

