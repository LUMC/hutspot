#!/usr/bin/env python3

import argparse
import sys


def main(args):
    bases = 0
    for reads, line in enumerate(sys.stdin, 1):
        bases += len(line.strip())
    with open(args.reads, 'w') as fout:
        print(reads, file=fout)
    with open(args.bases, 'w') as fout:
        print(bases, file=fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--reads', required=True)
    parser.add_argument('--bases', required=True)
    arguments = parser.parse_args()

    main(arguments)
