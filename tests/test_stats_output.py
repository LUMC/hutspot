#!/usr/bin/env python

import json
import pathlib
import pytest
import os


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_exists(workflow_dir):
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    assert os.path.exists(full_path)


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_name(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['sample_name'] == 'micro'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_mapped_reads(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['mapped_reads'] == '15424'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_mapped_bases(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['mapped_bases'] == '2262944'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_preqc_reads(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['preqc_reads'] == '15440'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_preqc_bases(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['preqc_bases'] == '2276743'


@pytest.mark.parametrize("json_field", [
    ('sample_name'), ('gender'), ('coverage'), ('picard_insertSize'),
    ('picard_AlignmentSummaryMetrics')])
@pytest.mark.workflow('test-integration-gene-bedfile')
def test_stats_file_preqc_bases(workflow_dir, json_field):
    """ Read in the stats json file """
    stats_file = 'stats.json'
    full_path = workflow_dir / pathlib.Path(stats_file)
    json_data = json.load(open(full_path))

    for sample in json_data['sample_stats']:
        assert json_field in sample
