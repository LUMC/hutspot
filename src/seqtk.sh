#!/usr/bin/env bash
set -e
set -o pipefail

count_json=${1}
input_fastq=${2}
output_fastq=${3}
max_bases=${4}

if [[ $max_bases -eq 'None' ]]; then
    ln -s $input_fastq $output_fastq
    exit 0
fi

if [[ -z $max_bases ]]; then
    ln -s $input_fastq $output_fastq
    exit 0
fi

bases=$(jq '.bases' $count_json)
frac=$(jq -n "$max_bases / $bases" | sed -e "s:e:E:g")
echo $frac
if (( $(echo "$frac > 1" | bc -l) )); then
    ln -s $input_fastq $output_fastq
else
    seqtk sample -s100 $frac $input_fastq | gzip -c > $output_fastq
fi
