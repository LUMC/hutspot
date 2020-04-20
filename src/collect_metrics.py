#!/usr/bin/env python3

import argparse
import json

from collections import defaultdict

def collect_cutadapt_summary(data, files):
    """ Read the different cutadapt files, sum the results and add the relevant
    metrics to the data dictionary """
    for filename in files:
        cutadapt = read_cutadapt(filename)
        data['preqc_reads'] += cutadapt['in_reads']
        data['preqc_bases'] += cutadapt['in_bp']
        data['postqc_reads'] += cutadapt['out_reads']
        # For some reason, cutadapt outputs the basepairs out separate for
        # forward and reverse reads
        data['postqc_bases'] += cutadapt['out_bp'] + cutadapt['out2_bp']

def read_cutadapt(filename):
    """ Read the cutadapt file and return the data """
    with open(filename, 'rt') as fin:
        header = next(fin).strip().split()
        values = next(fin).strip().split()

    data = dict()
    for field, value in zip(header, values):
        try:
            data[field] = int(value)
        except ValueError:
            data[field] = value

    return data

def main(args):
    # Data structure to store all data
    data = defaultdict(int)

    # Set the sample name
    data['sample_name'] = args.sample

    # Collect cutadapt data
    collect_cutadapt_summary(data, args.cutadapt_summary)

    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--cutadapt-summary',
                        required = True,
                        help = 'Cutadapt summary output',
                        nargs = '+'
    )
    parser.add_argument('--sample',
                        required = True,
                        help = 'Name of the sample'
    )

    args = parser.parse_args()

    main(args)
