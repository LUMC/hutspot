#!/usr/bin/env bash
set -eu
set -o pipefail

input_r1=${1}
input_r2=${2}
output_r1=${3}
output_r2=${4}
odir=${5}

fastqc --nogroup -o ${odir} ${input_r1} ${input_r2}

if [[ -f ${output_r1} ]]; then
    unzip -t ${output_r1} || truncate -s0 ${output_r1}
else
    touch ${output_r1}
fi

if [[ -f ${output_r2} ]]; then
    unzip -t ${output_r2} || truncate -s0 ${output_r2}
else
    touch ${output_r2}
fi
