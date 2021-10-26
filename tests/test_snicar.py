from funcs import generate_snicar_params_single_layer, call_snicar
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest


###########################
# SET USER DEFINE VARIABLES
###########################


savepath = '/home/joe/Code/BioSNICARParameterization/'
path_to_data = str(savepath+'snicar_data_single_layer.csv')
pytest.global_variable_1 = 0

ds = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/tests/matlab_benchmark_data.csv',header=None)

# @pytest.mark.parametrize("layer_type", [0,1])
# @pytest.mark.parametrize("density", [400, 500, 600, 700, 800])
# @pytest.mark.parametrize('reff', [200, 400, 600, 80, 1000])
# @pytest.mark.parametrize("zen", [30, 40, 50, 60, 70])
# @pytest.mark.parametrize("bc",[500, 1000, 1500, 2000])
# @pytest.mark.parametrize("dz",[[0.02,0.04,0.06, 0.08, 0.1],
#                                [0.04,0.06,0.08, 0.10, 0.15],
#                                [0.05,0.10,0.15, 0.2, 0.5],
#                                [0.15, 0.2, 0.25, 0.3, 0.5],
#                                [0.5,0.5,0.5,1,10]])#,
# def test_snicar_5_layer(layer_type, density, dz, zen, reff, bc):
#     params = generate_snicar_params_single_layer(layer_type, density, dz, 0, zen, reff, bc)             
#     albedo, BBA = call_snicar(params)
#     assert len(albedo)==480
#     assert sum(albedo)>0
#     assert BBA>0
#     assert abs(ds.iloc[480,pytest.global_variable_1] - BBA) < 0.002
#     pytest.global_variable_1+=1
    
#     ## READ IN MULTI LAYER CSV HERE, COMPARE ALBEDO COLUMN_WISE, PASS TEST IF DIFF<THRESHOLD


lyrList = [0,1]
densList = [400, 500, 600, 700, 800]
reffList = [200, 400, 600, 800, 1000]
zenList = [30, 40, 50, 60, 70]
bcList = [500, 1000, 1500, 2000]
dzList =[[0.02,0.04,0.06, 0.08, 0.1],
        [0.04,0.06,0.08, 0.10, 0.15],
        [0.05,0.10,0.15, 0.2, 0.5],
        [0.15, 0.2, 0.25, 0.3, 0.5],
        [0.5,0.5,0.5,1,10]]

specOut = np.zeros(shape=(5000,481))
errorList = []
counter = 0

for layer_type in lyrList:
    for density in densList:
        for reff in reffList:
            for zen in zenList:
                for bc in bcList:
                    for dz in dzList:

                        params = generate_snicar_params_single_layer(layer_type, density, dz, 0, zen, reff, bc)             
                        albedo, BBA = call_snicar(params)
                    
                        specOut[counter,0:480]  = albedo
                        specOut[counter,480] = BBA
                    
                        counter+=1

np.savetxt("test_results.csv", specOut, delimiter=",")


import numpy as np
import pandas as pd
py = pd.read_csv('test_results.csv', header=None)


bb_py = py.loc[:,481]
bb_mat = ds.loc[:,481]
error[error['abs_error']>1e-8].count()