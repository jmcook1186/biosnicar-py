from funcs import generate_snicar_params_single_layer, call_snicar
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest


###########################
# SET USER DEFINE VARIABLES
###########################

dzs = [0.1, 0.15, 0.2, 0.3, 0.5, 0.6, 0.7]
bc_conc = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
densities = [400, 500, 600, 700, 800, 850]
zeniths = [30, 40, 50, 60, 70, 80]
savepath = '/home/joe/Code/BioSNICARParameterization/'
path_to_data = str(savepath+'snicar_data_single_layer.csv')


@pytest.mark.parametrize("density", [400, 500])
@pytest.mark.parametrize("zen", [35, 45])
def test_snicar(density, zen):
    params = generate_snicar_params_single_layer(density, 1, 0, zen)             
    albedo, BBA = call_snicar(params)
    assert len(albedo)==480
    assert sum(albedo)>0
    assert BBA>0



