#!bin/bash

ANNOT_DIR="/Users/emantaha/Desktop/thesis/annotation"
CCRE_BED="/Users/emantaha/Desktop/thesis/cCRE/GRCh38-cCREs.bed"
RESULTS_DIR="/Users/emantaha/Desktop/thesis/results"
GTF="$ANNOT_DIR/gencode.v40.annotation.gtf"

# Extract all transcript TSSs from full v40

gawk '
BEGIN { OFS="\t" }

$3 == "transcript" {

    match($0, /gene_id "([^"]+)"/, gid)
    match($0, /gene_name "([^"]+)"/, gname)
    match($0, /gene_type "([^"]+)"/, gtype)

    tss = ($7 == "+") ? $4 : $5

    print $1, tss-1, tss, gid[1], gname[1], $7, gtype[1]
}
' "$GTF" | sort -k1,1 -k2,2n > "$ANNOT_DIR/TSS.Full.bed"


# create 2kb extended window

awk '{s=($2-2000<0)?0:$2-2000; print $1"\t"s"\t"$3+2000"\t"$4"\t"$6"\t"$7 }' \
    $ANNOT_DIR/TSS.FULL.bed | sort -k1,1 -k2,2n > $ANNOT_DIR/TSS_extended_2kb.Full.bed


# Extract PLS cCREs

awk '$6 == "PLS"' $CCRE_BED | sort -k1,1 -k2,2n > $RESULTS_DIR/PLS_cCREs_sorted.bed

# variables

REGIONS=$RESULTS_DIR/PLS_cCREs_sorted.bed
TSS=$ANNOT_DIR/TSS.Full.bed
PROX=$ANNOT_DIR/TSS_extended_2kb_Full.bed

# Intersect

bedtools intersect -wo -a $REGIONS -b $TSS > $RESULTS_DIR/tmp.tss
bedtools intersect -v -a $REGIONS -b $TSS > $RESULTS_DIR/tmp.no
bedtools intersect -u -a $RESULTS_DIR/tmp.no -b $PROX | \ 
    bedtools closet -d -a stdin -b $TSS > $RESULTS_DIR/tmp.distance



echo -e "cCRE_id\tgene_id\tgene_name\tgene_type\tdistance_bp" > $RESULTS_DIR/PLS_gene_assignment_full.txt

# distance

{
  awk 'BEGIN{OFS="\t"} {print $4, $10, $11, $13, 0}' $RESULTS_DIR/tmp.tss
  awk 'BEGIN{OFS="\t"} { m=int(($2+$3)/2); d=m-$9; if(d<0)d=-d; print $4, $10, $11, $13, d }' $RESULTS_DIR/tmp.distance
} | sort -u >> $RESULTS_DIR/PLS_gene_assignment_full.txt



 
