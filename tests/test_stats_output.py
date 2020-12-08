#!/usr/bin/env python

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

    assert data['mapped_reads'] == '15515'


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

    assert data['mapped_bases'] == '2275114'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_usable_reads(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['usable_reads'] == '15477'


@pytest.mark.workflow('test-integration-no-cluster')
def test_stats_file_usable_bases(workflow_dir):
    """ Read in the stats file and do some tests """
    stats_file = 'stats.tsv'
    full_path = workflow_dir / pathlib.Path(stats_file)
    # Since we only have two lines in the test stats file, this works
    with open(full_path, 'r') as fin:
        header = next(fin).strip().split('\t')
        values = next(fin).strip().split('\t')
    data = dict(zip(header, values))

    assert data['usable_bases'] == '2270739'
