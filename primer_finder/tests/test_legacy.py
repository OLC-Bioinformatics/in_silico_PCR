#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from spadespipeline.legacy_vtyper import Custom, Filer, Vtyper
from validator_helper import validate
from argparse import ArgumentParser
from time import time
import pytest
import shutil
import os

# Define global variables
test_path = os.path.abspath(os.path.dirname(__file__))
test_file_root = os.path.join(test_path, 'test_data')
test_sequences_path = os.path.join(test_file_root, 'sequences', 'legacy')
test_desired_output = os.path.join(test_file_root, 'desired_output')
test_primer_folder = os.path.join(test_file_root, 'primer_files')

__author__ = 'adamkoziol'


def create_args():
    arguments = ArgumentParser()
    arguments.sequencepath = test_sequences_path
    arguments.starttime = time()
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    arguments.runmetadata.samples = Filer.filer(arguments)
    return arguments


def create_validator(reference_csv, test_csv, name_column='Sample'):
    global validator, column_list
    column_list = validate.find_all_columns(csv_file=reference_csv,
                                            columns_to_exclude=[name_column])
    validator = validate.Validator(reference_csv=reference_csv,
                                   test_csv=test_csv,
                                   column_list=column_list,
                                   identifying_column=name_column)


def headers_same():
    # Check that the reference and test CSV files have the same headers.
    headers_are_same = validator.same_columns_in_ref_and_test()
    assert headers_are_same


def columns_present():
    # Check that all the columns in the column_list are in both the reference and test CSV files.
    columns = validator.all_test_columns_in_ref_and_test()
    assert columns


def samples_present():
    # Check that both test and reference have all the same sample names.
    samples = validator.check_samples_present()
    assert samples


def values():
    # Finally, check that Categorical columns are identical for test and reference, and range values
    # are within acceptable limits
    values_are_ok = validator.check_columns_match()
    assert values_are_ok


def clear_analyses():
    for root, dirs, files in os.walk(test_sequences_path):
        for directory in dirs:
            if os.path.isdir(os.path.join(root, directory)):
                shutil.rmtree(os.path.join(root, directory))


def remove_formatted_primers():
    shutil.rmtree(os.path.join(test_primer_folder, 'epcr_formatted_primers'))


def test_empty():
    with pytest.raises(TypeError):
        vt = Vtyper()


def test_custom():
    global args
    args = create_args()
    epcr = Custom(inputobject=args,
                  analysistype='custom_epcr',
                  primerfile=os.path.join(test_primer_folder, 'O157_LMhlyA_stn.txt'),
                  ampliconsize=10000,
                  mismatches=1,
                  export_amplicons=False)
    epcr.main()


def test_validate_custom():
    create_validator(reference_csv=os.path.join(test_desired_output, 'custom_legacy.csv'),
                     test_csv=os.path.join(test_sequences_path, 'reports', 'ePCR_report.csv'))
    headers_same()
    columns_present()
    samples_present()
    values()
    clear_analyses()
    remove_formatted_primers()


def test_vtyper():
    # Perform vtx typing
    vtyper = Vtyper(inputobject=args,
                    analysistype='vtyper_legacy',
                    mismatches=1)
    vtyper.vtyper()


def test_validate_vtyper():
    create_validator(reference_csv=os.path.join(test_desired_output, 'vtyper_legacy.csv'),
                     test_csv=os.path.join(test_sequences_path, 'reports', 'vtyper_legacy.csv'),
                     name_column='Strain')
    headers_same()
    columns_present()
    samples_present()
    values()
    clear_analyses()
