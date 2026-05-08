#!/bin/bash
#SBATCH --job-name=MEI_tissues_GAT
#SBATCH -p cbr_q_huge
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 5-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err


module load bedtools-2.30
module load bcftools-1.21
module load parallel
module load R

source ~/.bashrc
conda activate gat_env_py38

which gat-run.py >/dev/null || { echo "ERROR: gat-run.py not found"; exit 1; }

VCFS=(
  "/gpfs/data/user/shreyags/TE_work/filtered_vcfs/LINE1.filtered.vcf.gz"
  "/gpfs/data/user/shreyags/TE_work/filtered_vcfs/ALU.filtered.vcf.gz"
)

ROADMAP_DIR="/gpfs/data/user/shreyags/TE_work/data/roadmap_anno"
GENOME="/gpfs/data/user/shreyags/TE_work/filtered_vcfs/genome.bed"
OUT_DIR="/gpfs/data/user/shreyags/TE_work/replication_results/annotations/roadmap/enrichment_global_100k"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

ITERS=100000

build_freq_beds() {
  local VCF="$1" NAME="$2"

  bcftools view -f PASS "$VCF" \
    | bcftools +fill-tags -- -t AF \
    | bcftools query -f '%CHROM\t%POS\t%INFO/AF\n' \
    | awk '$1 ~ /^chr/ {print $1"\t"$2-1"\t"$2"\t"$3}' \
    > "${NAME}_all_with_af.bed"

  awk '{print $1"\t"$2"\t"$3}'               "${NAME}_all_with_af.bed" > "${NAME}.all.bed"
  awk '$4 > 0.01 {print $1"\t"$2"\t"$3}'     "${NAME}_all_with_af.bed" > "${NAME}.common.bed"

  rm -f "${NAME}_all_with_af.bed"
}

process_tissue_freq() {
  local ID="$1" FREQCLASS="$2" NAME="$3"

  local SEGFILE="${OUT_DIR}/${NAME}.${FREQCLASS}.bed"
  local ANNO_FILE="${ROADMAP_DIR}/${ID}_15_coreMarks_hg38lift_dense.bed.gz"
  local GTMPDIR="${OUT_DIR}/tmp_${NAME}_${ID}_${FREQCLASS}"

  [[ ! -s "$SEGFILE" ]] && return
  [[ ! -f "$ANNO_FILE" ]] && return

  mkdir -p "$GTMPDIR"

  local ANNO_BED="$GTMPDIR/anno.bed"

  # annotations building
  gzip -dc "$ANNO_FILE" \
    | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
    | sort -k1,1 -k2,2n \
    | awk '
      {
        key=$4
        if (NR==1 || key!=prev_key || $1!=prev_chr || $2>prev_end) {
          if (NR>1) print prev_chr"\t"prev_start"\t"prev_end"\t"prev_key
          prev_chr=$1; prev_start=$2; prev_end=$3; prev_key=key
        } else {
          if ($3>prev_end) prev_end=$3
        }
      }
      END {
        if (NR>0) print prev_chr"\t"prev_start"\t"prev_end"\t"prev_key
      }
    ' > "$ANNO_BED"

  # skip if annotation empty
  [[ ! -s "$ANNO_BED" ]] && return

  gat-run.py \
    --segments="$SEGFILE" \
    --annotations="$ANNO_BED" \
    --workspace="$GENOME" \
    --num-samples=$ITERS \
    --ignore-segment-tracks \
    > "$GTMPDIR/out.tsv" \
    2> "$GTMPDIR/gat.err"

  [[ ! -s "$GTMPDIR/out.tsv" ]] && return

  {
    flock 200
    awk -v tissue="$ID" \
      'NR>1 {print tissue"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$10"\t"$11}' \
      "$GTMPDIR/out.tsv" \
      >> "${OUT_DIR}/${NAME}_${FREQCLASS}.tsv"
  } 200>"${OUT_DIR}/${NAME}_${FREQCLASS}.lock"
}

export -f process_tissue_freq
export ROADMAP_DIR GENOME OUT_DIR ITERS

find "$ROADMAP_DIR" -name "*_15_coreMarks_hg38lift_dense.bed.gz" \
  | sed 's#.*/##;s/_15_coreMarks_hg38lift_dense.bed.gz//' \
  | sort > tissues.list

for VCF in "${VCFS[@]}"; do

  NAME=$(basename "$VCF" .vcf.gz)

  build_freq_beds "$VCF" "$NAME"

  for fc in all common; do
    echo -e "Tissue\tCategory\tObserved\tExpected\tFold\tlog2Fold\tPvalue\tQvalue_GAT" \
      > "${OUT_DIR}/${NAME}_${fc}.tsv"
  done

  parallel -j 12 \
    process_tissue_freq {1} {2} "$NAME" \
    :::: tissues.list \
    ::: all common

  rm -f ${OUT_DIR}/${NAME}_*.lock
  rm -rf ${OUT_DIR}/tmp_${NAME}_*

done

ls ${OUT_DIR}/*.tsv