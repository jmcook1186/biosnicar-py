#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr

"""Contains functions for updating optical property files and formatting for reading in to BioSNICAR.
"""

ds = xr.open_dataset(
    "/Data/GO_files/dust_greenland_C_20150308.nc"
)

wavelengths = np.arange(0.2, 5, 0.01)

# initial values
ext_cff_mss = list(ds.ext_cff_mss.values)
ss_alb = list(ds.ss_alb.values)
asm_prm = list(ds.asm_prm.values)

# set up output array
out = np.zeros(shape=(len(wavelengths), 3))

# initialise counter
COUNT = 0

# loop through relevant vars
# and use the first value to
# pad 10 elements at front of array
for y in [ext_cff_mss, ss_alb, asm_prm]:
    num_elements_to_add = 10
    pad_value = y[0]
    y = [pad_value for _ in range(10)] + y
    out[:, COUNT] = y
    COUNT += 1

# alternative = use interp1d
# x = ds.wvl.values
# for y in [ext_cff_mss, ss_alb, asm_prm]:

#     xnew = wavelengths*1e-6
#     f = interp1d(x,y,kind='cubic',fill_value=y[0])
#     ynew = f(xnew)
#     out[:,count] = ynew
#     count +=1

# update variables with newly extended values
ext_cff_mss = out[:, 0]
ss_alb = out[:, 1]
asm_prm = out[:, 2]

# create new xarray dataset
newDS = xr.Dataset(
    data_vars=dict(
        ext_cff_mss=(["wvl"], ext_cff_mss),
        ss_alb=(["wvl"], ss_alb),
        asm_prm=(["wvl"], asm_prm),
    ),
    coords=dict(wvl=(["wvl"], wavelengths)),
)

# save new netCDF file
newDS.to_netcdf(
    "/home/joe/Code/BioSNICAR_GO_PY/Data/GO_files/480band/Cook_Greenland_Central.nc",
    mode="w",
)
