#!/bin/bash

RESULTS_DIR="/Users/emantaha/Desktop/thesis/results"
BEDS_DIR="$RESULTS_DIR/beds"
OUT_DIR="$RESULTS_DIR/geodist_direction"

mkdir -p "$OUT_DIR"

for mark in H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K9me3; do

for pos_bed in "$BEDS_DIR"/pos_"${mark}"_*.bed; do

```
[[ -f "$pos_bed" ]] || continue

met=$(basename "$pos_bed")
met=${met#pos_${mark}_}
met=${met%.bed}

for group in pos neg; do

  bed="$BEDS_DIR/${group}_${mark}_${met}.bed"
  out="$OUT_DIR/dist_${group}_${mark}_${met}.txt"

  [[ -f "$bed" ]] || continue
  [[ $(wc -l < "$bed") -lt 2 ]] && continue

  sort -k1,1 -k2,2n "$bed" |
  bedtools closest -a - -b "$bed" -d -io |
  awk 'BEGIN{OFS="\t"}
       NR==1{print "cCRE_1","cCRE_2","chr","distance_bp"}
       $NF>=0 && $8!="."{print $4,$8,$1,$NF}' \
  > "$out"

done
```

done
done

echo "Done"
