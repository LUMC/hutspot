#!/usr/bin/env bash

count_json=${1}
input_fastq=${2}
output_fastq=${3}
max_bases=${4}


bases=$(jq '.bases' $count_json)
frac=$(jq -n "$max_bases / $bases" | sed -e "s:e:E:g")
echo $frac
if (( $(echo "$frac > 1" | bc -l) )); then
    ln -s $input_fastq $output_fastq
else
    seqtk sample -s100 $frac $input_fastq | gzip -c > $output_fastq
fi
