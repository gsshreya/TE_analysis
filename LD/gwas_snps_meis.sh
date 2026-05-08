#!/bin/bash
#SBATCH --job-name=gwas_snps_meis
#SBATCH -p cbr_q_small
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 05:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

cd $HOME/TE_work/replication_results/ld

GWAS="gwas_chrpos.txt"
LD="ld_mei_snp.ld"
OUTDIR="$HOME/TE_work/replication_results/ld/GWAS"
mkdir -p $OUTDIR

for TYPE in ALU LINE1 SVA
do
    SNPFILE="${TYPE}_snps.txt"

    echo "Processing $TYPE"

    # Find SNPs overlapping GWAS loci
    grep -xFf <(cut -d':' -f1,2 "$SNPFILE") "$GWAS" \
        | sort -u > $OUTDIR/${TYPE}_GWAS_snps.txt

    COUNT=$(wc -l < $OUTDIR/${TYPE}_GWAS_snps.txt)

    echo "$TYPE GWAS SNP overlap: $COUNT" | tee -a $OUTDIR/GWAS_overlap_counts.txt

    # Extract chr:pos:ref:alt versions of those SNPs
    grep -Ff $OUTDIR/${TYPE}_GWAS_snps.txt "$SNPFILE" \
        > $OUTDIR/${TYPE}_GWAS_snps_fullID.txt

    # Find which MEIs those SNPs are in LD with
    awk -v type="$TYPE" '
    NR==FNR {gwas[$0]; next}

    {
        snp1=$3
        snp2=$6
        mei=""

        if($3 ~ /INS:ME:/) mei=$3
        if($6 ~ /INS:ME:/) mei=$6

        pos1=substr(snp1,1,index(snp1,":")-1)
        pos2=substr(snp2,1,index(snp2,":")-1)

        split(snp1,a,":")
        split(snp2,b,":")

        chrpos1=a[1]":"a[2]
        chrpos2=b[1]":"b[2]

        if(chrpos1 in gwas && mei ~ type) print chrpos1,mei
        if(chrpos2 in gwas && mei ~ type) print chrpos2,mei
    }
    ' $OUTDIR/${TYPE}_GWAS_snps.txt "$LD" \
    | sort -u > $OUTDIR/${TYPE}_GWAS_MEI_pairs.txt

done
