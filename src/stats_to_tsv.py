import argparse
import json

from collections import OrderedDict
from pathlib import Path


def get_vcf_stats(sample_name, vcfstats):
    vcf_sample = next(x for x in vcfstats['samples'] if x['name'] == sample_name)
    return {
        "total_variants": vcf_sample['total_variants'],
        "snps": vcf_sample['variant_types']['snps'],
        "insertions": vcf_sample['variant_types']['insertions'],
        "deletions": vcf_sample['variant_types']['deletions'],
        "transversions": vcf_sample['transversions'],
        "transitions": vcf_sample['transitions'],
        "ti_tv_ratio": vcf_sample['ti_tv_ratio'],
        "homozygous_variants": vcf_sample['genotypes']['hom_alt'],
        "heterozygous_variants": vcf_sample['genotypes']['het']
    }


def get_covstats(cov_d):
    s_d = cov_d['covstats']['stats']['coverage']['_all']
    tmp_d = {
        'median_coverage': s_d['median'],
        'mean_coverage': s_d['mean'],
        'modal_coverage': s_d['mode'],
        'horizontal_coverage': s_d['horizontal'],
        'coverage_frac_min_100x': s_d['frac_min_100x'],
        'coverage_frac_min_10x': s_d['frac_min_10x'],
        'coverage_frac_min_20x': s_d['frac_min_20x'],
        'coverage_frac_min_30x': s_d['frac_min_30x'],
        'coverage_frac_min_40x': s_d['frac_min_40x'],
        'coverage_frac_min_50x': s_d['frac_min_50x'],
        'determined_gender': cov_d['gender']
    }
    return {"{0}_{1}".format(cov_d['name'], k): v for k, v in tmp_d.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to input JSON", type=Path)

    args = parser.parse_args()

    with args.input.open() as handle:
        orig_dict = json.load(handle)

    sdicts = []

    vcfstats = orig_dict['multisample_vcfstats']

    for sample in orig_dict['sample_stats']:
        sname = sample['sample_name']
        sample_dict = OrderedDict()
        sample_dict.update({
            "sample_name": sname,
            "preqc_reads": sample['pre_qc_fastq_count']['reads'],
            "preqc_bases": sample['pre_qc_fastq_count']['bases'],
            "postqc_reads": sample['post_qc_fastq_count']['reads'],
            "postqc_bases": sample['post_qc_fastq_count']['bases'],
            "mapped_reads": sample['n_mapped_reads'],
            "mapped_bases": sample['n_mapped_bases'],
            "usable_reads": sample['n_usable_reads'],
            "usable_bases": sample['n_usable_bases']
        })
        sample_dict.update(get_vcf_stats(sname, vcfstats))
        if "covstats" in sample:
            for cov_d in sample['covstats']:
                sample_dict.update(get_covstats(cov_d))
        sdicts.append(sample_dict)

    lens = [len(list(x.keys())) for x in sdicts]
    longest_dict = sdicts[lens.index(max(lens))]

    header = "\t".join(longest_dict.keys())
    print(header)
    for sd in sdicts:
        vals = []
        for k in longest_dict.keys():
            if k in sd:
                vals.append(sd[k])
            else:
                vals.append("NA")
        print("\t".join(map(str, vals)))