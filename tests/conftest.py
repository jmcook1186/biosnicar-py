import pandas as pd
import pytest


@pytest.fixture
def get_matlab_data():
    return pd.read_csv("./tests/test_data/matlab_benchmark_data.csv", header=None)


@pytest.fixture
def get_matlab_data_toon():
    return pd.read_csv("./tests/test_data/matlab_benchmark_data_toon.csv", header=None)


@pytest.fixture
def get_python_data():
    return pd.read_csv("./tests/test_data/py_benchmark_data.csv", header=None).transpose()


@pytest.fixture
def get_python_data_toon():
    return pd.read_csv("./tests/test_data/py_benchmark_data_toon.csv", header=None).transpose()


@pytest.fixture
def set_tolerance():
    return 1e-6


@pytest.fixture
def get_n_spectra():
    return 25


@pytest.fixture
def fuzz():
    return True

@pytest.fixture
def new_benchmark_ad():
    return True

@pytest.fixture
def new_benchmark_toon():
    return True