import click
import cyvcf2
import numpy
import json

from tqdm import tqdm

from collections import Counter


def gen_chrom_counter(vcf_reader):
    """Generate chromosome counter from VCF reader"""
    return Counter({n: 0 for n in vcf_reader.seqnames})


class Sample(object):
    def __init__(self, name, idx):
        self.name = name
        self.idx = idx
        self.transversions = 0
        self.transitions = 0
        self.hom_ref = 0
        self.het = 0
        self.hom_alt = 0
        self.deletions = 0
        self.insertions = 0
        self.snps = 0

        self.__gq_counter = Counter({i: 0 for i in range(100)})

    @property
    def ti_tv(self):
        if self.transversions > 0:
            return float(self.transitions)/self.transversions
        return numpy.nan

    def add_variant(self, var):
        """assuming gts012=True"""
        typ = var.gt_types[self.idx]
        if typ == 3:
            return None

        gq = var.gt_quals[self.idx]
        self.__gq_counter[gq] += 1

        if typ == 0:
            self.hom_ref += 1
            return None

        if typ == 1:
            self.het += 1
        if typ == 2:
            self.hom_alt += 1

        if var.is_snp and var.is_transition:  # this only works in python 2 for now. See: https://github.com/brentp/cyvcf2/pull/70
            self.transitions += 1
            self.snps += 1
        elif var.is_snp:
            self.transversions += 1
            self.snps += 1
        elif var.is_indel and var.is_deletion:
            self.deletions += 1
        elif var.is_indel:
            self.insertions += 1

    @property
    def gq_distr(self):
        return [self.__gq_counter[x] for x in range(100)]

    @property
    def total_variants(self):
        return self.hom_alt + self.het

    @property
    def as_dict(self):
        return {
            "name": self.name,
            "total_variants": self.total_variants,
            "variant_types": {
                "snps": self.snps,
                "deletions": self.deletions,
                "insertions": self.insertions,
            },
            "genotypes": {
                "hom_ref": self.hom_ref,
                "het": self.het,
                "hom_alt": self.hom_alt
            },
            "transitions": self.transitions,
            "transversions": self.transversions,
            "ti_tv_ratio": self.ti_tv,
            "gq_distribution": self.gq_distr
        }


class Stats(object):
    def __init__(self, vcf_path):
        self.path = vcf_path
        self.vcf = cyvcf2.VCF(vcf_path, gts012=True)

        self.samples = [Sample(x, i) for i, x in enumerate(self.vcf.samples)]
        self.chrom_counter = gen_chrom_counter(self.vcf)

        self.__calculated = False

    def calculate(self):
        for record in tqdm(self.vcf, unit="variants", unit_scale=True):
            for s in self.samples:
                s.add_variant(record)
            self.chrom_counter[record.CHROM] += 1
        self.__calculated = True

    @property
    def total_variants(self):
        return sum(self.chrom_counter.values())

    @property
    def as_dict(self):
        if not self.__calculated:
            self.calculate()
        return {
            "vcf_path": self.path,
            "total_variants": self.total_variants,
            "samples": [s.as_dict for s in self.samples],
            "per_chromosome_variants": {k: v for k, v in self.chrom_counter.items()}
        }

    @property
    def as_json(self):
        return json.dumps(self.as_dict, sort_keys=True, indent=4)


@click.command()
@click.option("-i",
              "--input",
              type=click.Path(exists=True, dir_okay=False, readable=True),
              required=True,
              help="Input VCF file")
def main(input):
    stats = Stats(input)
    print(stats.as_json)


if __name__ == "__main__":
    main()
