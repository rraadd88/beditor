"""
Test runs.
Notes: The paths are relative to '../'
"""
import logging
import os
from beditor.run import cli

def set_wd():
    print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())
    if not os.getcwd().endswith("/examples"):
        os.chdir('./examples/')
        print(">>>>>>>>>>>>>>>>> current dir:"+os.getcwd())

def test_base_positions():
    set_wd()
    return cli(
        config_path='inputs/mutations/base/positions.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )

def test_base_regions():
    set_wd()
    return cli(
        config_path='inputs/mutations/base/regions.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )

def test_base_regions():
    set_wd()
    return cli(
        config_path='inputs/mutations/base/points.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )

def test_protein_positions():
    set_wd()
    return cli(
        config_path='inputs/mutations/protein/positions.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )

def test_protein_regions():
    set_wd()
    return cli(
        config_path='inputs/mutations/protein/regions.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )

def test_protein_regions():
    set_wd()
    return cli(
        config_path='inputs/mutations/protein/points.yml',
        verbose='INFO',#'DEBUG',
        skip=['Visualize'], # skipped because interactive
        force=True,
    )