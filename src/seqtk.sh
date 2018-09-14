#!/usr/bin/env bash
set -ex
set -o pipefail

count_json=${1}
input_fastq=${2}
output_fastq=${3}
max_bases=${4}

if [[ ${max_bases} -eq 'None' ]]; then
    ln -s ${input_fastq} ${output_fastq}
    exit 0
fi

if [[ -z $max_bases ]]; then
    ln -s ${input_fastq} ${output_fastq}
    exit 0
fi

bases=$(jq '.bases' $count_json)
frac=$(jq -n "$max_bases / $bases" | sed -e "s:e:E:g")
echo $frac
frac_higher_than_one=$(echo "${frac} > 1" | bc -l)
if [[ ${frac_higher_than_one} -eq 1 ]]; then
    ln -s ${input_fastq} ${output_fastq}
else
    seqtk sample -s100 ${input_fastq} ${frac} | gzip -c > ${output_fastq}
fi
