import pandas as pd
import pytest


@pytest.fixture
def get_matlab_data():
    return pd.read_csv("./tests/test_data/matlab_benchmark_data.csv", header=None)


@pytest.fixture
def get_matlab_data_clean():
    return pd.read_csv("./tests/test_data/matlab_benchmark_data_clean.csv", header=None)

@pytest.fixture
def get_matlab_data_toon():
    return pd.read_csv("./tests/test_data/matlab_benchmark_data_toon.csv", header=None)


@pytest.fixture
def get_python_data():
    return pd.read_csv("./tests/test_data/py_benchmark_data.csv", header=None).transpose()


@pytest.fixture
def get_python_data_clean():
    return pd.read_csv("./tests/test_data/py_benchmark_data_clean.csv", header=None).transpose()


@pytest.fixture
def get_python_data_toon():
    return pd.read_csv("./tests/test_data/py_benchmark_data_toon.csv", header=None).transpose()


@pytest.fixture
def set_tolerance():
    return 1e-3


@pytest.fixture
def get_n_spectra():
    return 25


@pytest.fixture
def fuzz():
    return False

@pytest.fixture
def new_benchmark_ad():
    return False

@pytest.fixture
def new_benchmark_ad_clean():
    return False

@pytest.fixture
def new_benchmark_toon():
    return False
