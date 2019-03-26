#!/usr/bin/env bash
#   hutspot - a DNAseq variant calling pipeline
#   Copyright (C) 2017-2019, Sander Bollen, Leiden University Medical Center
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
set -ex
set -o pipefail

count_json=${1}
input_fastq=${2}
output_fastq=${3}
max_bases=${4}

full_input=$(readlink -f ${input_fastq})
if [[ ${max_bases} -eq 'None' ]]; then
    ln -s ${full_input} ${output_fastq}
    exit 0
fi

if [[ -z $max_bases ]]; then
    ln -s ${full_input} ${output_fastq}
    exit 0
fi

bases=$(jq '.bases' $count_json)
frac=$(jq -n "$max_bases / $bases" | sed -e "s:e:E:g")
echo $frac
frac_higher_than_one=$(echo "${frac} > 1" | bc )
if [[ ${frac_higher_than_one} -eq 1 ]]; then
    ln -s ${full_input} ${output_fastq}
else
    seqtk sample -s100 ${input_fastq} ${frac} | gzip -c > ${output_fastq}
fi
