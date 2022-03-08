#!/usr/bin/python

import pandas as pd
import pytest


@pytest.fixture
def get_matlab_data():
    """Loads benchmark data from matlab model.
    """
    return pd.read_csv("./tests/test_data/matlab_benchmark_data.csv", header=None)


@pytest.fixture
def get_matlab_data_clean():
    """Loads benchmark data from matlab model for clean ice.
    """
    return pd.read_csv("./tests/test_data/matlab_benchmark_data_clean.csv", header=None)


@pytest.fixture
def get_python_data():
    """Loads data from BioSNICAR for benchmarking.
    """
    return pd.read_csv(
        "./tests/test_data/py_benchmark_data.csv", header=None
    ).transpose()


@pytest.fixture
def get_python_data_clean():
    """Loads data from BioSNICAR for benchmarking clean ice.
    """
    return pd.read_csv(
        "./tests/test_data/py_benchmark_data_clean.csv", header=None
    ).transpose()


@pytest.fixture
def set_tolerance():
    """Sets the error tolerance for tests to pass. Default 1e-5
    """
    return 1e-5


@pytest.fixture
def get_n_spectra():
    """Defines how many random spectra pairs to plot.
    """
    return 25


@pytest.fixture
def fuzz():
    """Toggles fuzzing tests on/off.
    """
    return True


@pytest.fixture
def new_benchmark_ad():
    """Toggles generation of new BioSNICAR benchmarking data on/off.
    """
    return True


@pytest.fixture
def new_benchmark_ad_clean():
    """Toggles generation of new BioSNICAR becnhmarking data on/off for clean ice.
    """
    return False


if __name__ == '__main__':
    pass