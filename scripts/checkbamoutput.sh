#!/usr/bin/env sh
input="$1"
output="$2"
outdir="$3"

if grep -Fxq "No errors found" $input; then
    echo "BAM file for ${outdir} good" > $output
else
    echo "Bad BAM file for ${outdir}" > $output
    mv ${outdir}/bt2.sorted.bam ${outdir}/bad_bt2.sorted.bam
fi
