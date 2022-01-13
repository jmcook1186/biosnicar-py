from benchmarking_funcs import generate_snicar_params, call_snicar
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import random

"""
To run configure these tests, update the values in conftest.py
Then navigate to the tests folder and run

`pytest`

The tests will automatically run - green dots indicate tests 
passing successfully. A plot of N random spectra pairs will
be saved to the /tests folder.

"""

def test_realistic_BBA(get_matlab_data, get_python_data):
    # are the values predicted by the model always physical (i.e. between 0-1)
    # do the files have the right shape and size?

    mat = get_matlab_data
    py = get_python_data
    
    bb_py = py.loc[:,481]
    bb_mat = mat.loc[:,481]

    assert len(bb_py) == len(bb_mat)
    assert bb_py[bb_py>1].count() ==0 and bb_py[bb_py<0].count() ==0
    assert bb_mat[bb_mat>1].count() ==0 and bb_mat[bb_mat<0].count() ==0
    

def test_compare_pyBBA_to_matBBA(get_matlab_data, get_python_data, set_tolerance):
    # check the BBA predicted for each run matches to within tolerance between
    # the two models
        mat = get_matlab_data
        py = get_python_data
        tol = set_tolerance
        bb_py = py.loc[:,481]
        bb_mat = mat.loc[:,481]
        error = np.array(abs(bb_mat -bb_py))
        assert  len(error[error>tol]) ==0


def test_compare_pySPEC_to_matSPEC(get_matlab_data, get_python_data, set_tolerance):
    # check that for each individual wavelenght, the spectral albedo
    # matches to within tolerance between the two models
        mat = get_matlab_data
        py = get_python_data
        tol = set_tolerance
        bb_spec = py.loc[:,:480]
        bb_spec = mat.loc[:,:480]
        error = np.array(abs(bb_spec -bb_spec))
        assert  len(error[error>1e-8]) ==0


def test_plot_random_spectra_pairs(get_matlab_data, get_python_data, get_n_spectra):
        # grabs n spectra and plots the python and matlab versions
        # for visual comparison
        mat = get_matlab_data
        py = get_python_data
        n_spec = get_n_spectra
        idxs = random.sample(range(0, py.shape[0]), n_spec)

        wavelength = np.arange(200,5000,10)
        py_spectra = py.iloc[0:-1,idxs]
        mat_spectra = mat.iloc[0:-1,idxs]

        plt.plot(wavelength, py_spectra) 
        plt.plot(wavelength, mat_spectra, linestyle=None,\
            marker='x')
        plt.xlabel("wavelength (nm)")
        plt.ylabel('Albedo')
        plt.title("solid lines = Python\ncrosses = Matlab")
        
        plt.savefig('py_mat_comparison.png')


