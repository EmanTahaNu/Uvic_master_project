#!/bin/bash

# Paths
SSD_ROOT="/Volumes/eman/thesis"
BIGWIG_DIR="$SSD_ROOT/bigwigs"
CCRE_BED="$SSD_ROOT/GRCh38-cCREs.bed"
CCRE_BED4="$SSD_ROOT/GRCh38-cCREs_4col.bed"
OUT_DIR="$SSD_ROOT/results/chipseq_signal"
TOOL="$SSD_ROOT/bigWigAverageOverBed"
# Create output folder
mkdir -p "$OUT_DIR"

# Create 4-column BED
awk '{print $1"\t"$2"\t"$3"\t"NR}' "$CCRE_BED" > "$CCRE_BED4"


# Run for every bigWig file
for bw in "$BIGWIG_DIR"/*.bigWig; do
    sample=$(basename "$bw" .bigWig)
    output="$OUT_DIR/${sample}_signal.tab"


    "$TOOL" "$bw" "$CCRE_BED4" "$output"

    if [ -f "$output" ]; then
        echo "Finished $sample"
    else
        echo "Failed $sample"
    fi
done

