#!/usr/bin/env python3

import argparse
import json

def print_tsv(data, filename):
    """ Flatten data and print in tsv format """
    # Determine the header
    for sample in data:
        header = sorted(data[sample].keys())

    # Print the header, and then for each sample the data
    with open(filename, 'w') as fout:
        print('sample_name', *header, sep='\t', file=fout)
        for sample in data:
            print(sample, *(data[sample][field] for field in header), sep='\t',
                 file=fout)


def main(args):
    data = dict()
    for filename in args.metrics:
        with open(filename, 'rt') as fin:
            metrics = json.load(fin)
            name = metrics.pop('sample_name')
            data[name] = metrics

    if args.json:
        with open(args.json, 'wt') as fout:
            json.dump(data, fout, indent=2)
    if args.tsv:
        print_tsv(data, args.tsv)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--metrics',
                        required = True,
                        help = 'Metrics json files',
                        nargs = '+'
    )
    parser.add_argument('--json', required=False)
    parser.add_argument('--tsv', required=False)

    args = parser.parse_args()

    main(args)
