#!/bin/bash
#SBATCH --job-name=te_annovar
#SBATCH -p cbr_q_small
#SBATCH --mem=32G
#SBATCH -c 32
#SBATCH -o logs/te_annovar_%j.out
#SBATCH -e logs/te_annovar_%j.err

module load bcftools-1.21


ANNOVAR_DB="/gpfs/data/user/shreyags/TE_work/annovar/humandb"
BUILD=hg38
VCF_DIR="$HOME/TE_work/filtered_vcfs/"
OUT_DIR="$HOME/TE_work/replication_results/insertions_annotations"
mkdir -p "$OUT_DIR"
SUMMARY_FILE="$OUT_DIR/annotation_summary.txt"
> "$SUMMARY_FILE"

for MEI in ALU LINE1; do

  echo "Converting to avinput"
  bcftools query \
  -f '%CHROM\t%POS\t%POS\t0\t-\n' \
  "$VCF_DIR/${MEI}.filtered.final.vcf.gz" > "$OUT_DIR/${MEI}.avinput"

  echo "Annotating ${MEI} MEIs..."

    
  /gpfs/data/user/shreyags/TE_work/annovar/table_annovar.pl \
    "$OUT_DIR/${MEI}.avinput" ${ANNOVAR_DB} \
    -buildver ${BUILD} \
    -out $OUT_DIR/${MEI}_ensGene_annovar \
    -protocol ensGene \
    -operation g \
    -nastring . \
    -polish \
    -arg '-splicing 5'
    
  echo "${MEI} - ensGene" >> "$SUMMARY_FILE"
    
    
  awk 'NR>1 {count[$6]++} END{for (c in count) print c, count[c]}' \
      "$OUT_DIR/${MEI}_ensGene_annovar.hg38_multianno.txt" >> "$SUMMARY_FILE"
      
  echo "" >> "$SUMMARY_FILE"
      
      
#  /gpfs/data/user/shreyags/TE_work/annovar/table_annovar.pl \
#    "$OUT_DIR/${MEI}.avinput" ${ANNOVAR_DB} \
#    -buildver ${BUILD} \
#    -out $OUT_DIR/${MEI}_refGeneWithVer_annovar \
#    -protocol refGeneWithVer \
#    -operation g \
#    -nastring . \
#    -polish \
#    -arg '-splicing 5'
    
#  echo "${MEI} - refGene" >> "$SUMMARY_FILE"
    
    
#  awk 'NR>1 {count[$6]++} END{for (c in count) print c, count[c]}' \
#      "$OUT_DIR/${MEI}_refGeneWithVer_annovar.hg38_multianno.txt" >> "$SUMMARY_FILE"
      
#  echo "" >> "$SUMMARY_FILE"


done

