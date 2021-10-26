import pandas as pd
import pytest

@pytest.fixture
def get_matlab_data():
    return pd.read_csv('/matlab_benchmark_data.csv',header=None)

@pytest.fixture
def get_python_data():
    return pd.read_csv('py_benchmark_data.csv', header=None).transpose()

@pytest.fixture
def set_tolerance():
    return 1e-8

@pytest.fixture
def get_n_spectra():
    return 25