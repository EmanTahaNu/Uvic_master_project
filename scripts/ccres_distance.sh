RESULTS_DIR="/Users/emantaha/Desktop/thesis/results"

# bedtool closest 

# -io ignores self overlaps
# -a, -b (for each region in -a ,find its nearst neighbor in -b)

bedtools closest \
-a $RESULTS_DIR/PLS_cCREs_sorted.bed\
-b $RESULTS_DIR/PLS_cCREs_sorted.bed\
-d -io > $RESULTS_DIR/tmp_pairs.txt

# extract needed columns

echo -e "cCRE_1\tcCRE_2\tchr\tdistance_bp" > $RESULTS_DIR/PLS_pairs.txt

awk 'BEGIN{OFS="\t"} {print $4, $10, $1, $NF}' \
    $RESULTS_DIR/tmp_pairs.txt >> $RESULTS_DIR/PLS_pairs.txt

 
