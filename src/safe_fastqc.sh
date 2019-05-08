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
set -eu
set -o pipefail

input_r1=${1}
input_r2=${2}
output_r1=${3}
output_r2=${4}
odir=${5}

fastqc --nogroup -o ${odir} ${input_r1} ${input_r2}

if [[ -f ${output_r1} ]]; then
    unzip -l ${output_r1} || truncate -s0 ${output_r1}
else
    touch ${output_r1}
fi

if [[ -f ${output_r2} ]]; then
    unzip -l ${output_r2} || truncate -s0 ${output_r2}
else
    touch ${output_r2}
fi
