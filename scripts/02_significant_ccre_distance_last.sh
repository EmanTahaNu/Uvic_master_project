#!/bin/bash


BASE="/Volumes/eman/thesis"
SIG_DIR="$BASE/results/FDR_0.1/significant"
OUT_DIR="$BASE/results"

echo " Extract significant cCRE IDs "
cat $SIG_DIR/PLS_correlation_*_significant.csv | \
    awk -F',' '{print $1}' | sort -u > $OUT_DIR/sig_ccre_ids_all.txt
echo "Unique significant cCRE IDs: $(wc -l < $OUT_DIR/sig_ccre_ids_all.txt)"

echo " Extract coordinates from cCRE reference"
awk 'NR==FNR{ids[$1]; next} $4 in ids' \
    $OUT_DIR/sig_ccre_ids_all.txt \
    $BASE/cCRE/GRCh38-cCREs.bed > $OUT_DIR/sig_ccre_coords.bed
echo "Coordinates extracted: $(wc -l < $OUT_DIR/sig_ccre_coords.bed)"



echo " Sort coordinates"
sort -k1,1 -k2,2n $OUT_DIR/sig_ccre_coords.bed > $OUT_DIR/sig_ccre_coords_sorted.bed



echo " Compute nearest-neighbor distances (bedtools closest) ==="
bedtools closest \
    -a $OUT_DIR/sig_ccre_coords_sorted.bed \
    -b $OUT_DIR/sig_ccre_coords_sorted.bed \
    -d -io > $OUT_DIR/tmp_sig_pairs.txt

echo -e "cCRE_1\tcCRE_2\tchr\tdistance_bp" > $OUT_DIR/sig_PLS_pairs.txt
awk 'BEGIN{OFS="\t"} {print $4, $10, $1, $NF}' \
    $OUT_DIR/tmp_sig_pairs.txt >> $OUT_DIR/sig_PLS_pairs.txt

echo "Pairs written: $(wc -l < $OUT_DIR/sig_PLS_pairs.txt)"

echo "Unique cCREs in pairs file: $(awk 'NR>1{print $1}' $OUT_DIR/sig_PLS_pairs.txt | sort -u | wc -l)"
echo "Unique cCREs in ID list:    $(wc -l < $OUT_DIR/sig_ccre_ids_all.txt)"
echo "Distance summary (bp):"
awk 'NR>1{print $NF}' $OUT_DIR/sig_PLS_pairs.txt | sort -n | \
    awk '{a[NR]=$1} END{print "  min:",a[1]; print "  median:",a[int(NR/2)]; print "  max:",a[NR]}'

