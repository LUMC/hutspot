import argparse
from pathlib import Path
import zipfile


class FastqcModule(object):
    """Helper class for Fastqc Modules"""

    def __init__(self, name):
        self.name = name
        self.header = None
        self.rows = []

    def add_header(self, headerstring):
        self.header = headerstring.split("\t")


def extract_data_txt(zip_path: Path, encoding='utf-8') -> str:
    """Extract text of data.txt from fastqc zip file"""
    with zipfile.ZipFile(zip_path) as zp:
        iflist = zp.infolist()
        data_info = next(x for x in iflist if x.filename.endswith("fastqc_data.txt"))
        with zp.open(data_info.filename) as data_handle:
            return data_handle.read().decode(encoding)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--preqc-r1",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for forward strand of merged fastq file before pre-processing")
    parser.add_argument("--preqc-r2",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for reverse strand of merged fastq file before pre-processing")
    parser.add_argument("--postqc-r1",
                        type=Path,
                        required=True,
                        help="Fastqc zip file for forward strand of fastq file after pre-processing")
    parser.add_argument("--postqc-r2",
                        type=Path,
                        required=True,
                        help="Fastqc zip for file reverse strand of fastq file after pre-processing")
    args = parser.parse_args()


