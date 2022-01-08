from benchmarking_funcs import generate_snicar_params, call_snicar
import numpy as np
import pandas as pd

lyrList = [0,1]
densList = [400, 500, 600, 700, 800]
reffList = [200, 400, 600, 800, 1000]
zenList = [30, 40, 50, 60, 70]
bcList = [0, 1000, 2000]
dust1List = [0, 10000, 20000, 50000]
dust5List = [0, 10000, 20000, 50000]
dzList =[[0.02,0.04,0.06, 0.08, 0.1],
        [0.04,0.06,0.08, 0.10, 0.15],
        [0.05,0.10,0.15, 0.2, 0.5],
        [0.15, 0.2, 0.25, 0.3, 0.5],
        [0.5,0.5,0.5,1,10]]

ncols = len(lyrList)*len(densList)*len(reffList)*len(zenList)\
    *len(bcList)*len(dzList)*len(dust1List)*len(dust5List)
print(ncols)
specOut = np.zeros(shape=(ncols,481))
errorList = []
counter = 0

for layer_type in lyrList:
    for density in densList:
        for reff in reffList:
            for zen in zenList:
                for bc in bcList:
                    for dz in dzList:
                        for dust1 in dust1List:
                            for dust5 in dust5List:

                                params = generate_snicar_params(layer_type, density, dz, 0, zen, reff, bc, dust1, dust5)             
                                albedo, BBA = call_snicar(params)
                            
                                specOut[counter,0:480]  = albedo
                                specOut[counter,480] = BBA
                            
                                counter+=1

np.savetxt("test_results.csv", specOut, delimiter=",")
