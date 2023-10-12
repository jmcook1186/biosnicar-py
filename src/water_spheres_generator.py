#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 22:38:08 2023

@author: Lou-Anne Chevrollier, University of Aarhus
"""

import pandas as pd
import numpy as np
from miepython import mie
import xarray as xr
import glob

def calc_mie_params(radius, n_in, k_in, n_ext):
    wvl = np.arange(0.205, 5, 0.01) # um
    X = 2 * np.pi * radius / wvl  # unitless
    n = n_in / n_ext
    k = k_in / n_ext
    qext, qsca, qback, g = mie(n - 1j * k, X) # cross sections in um2 
    qabs = qext - qsca
    assym = g
    ss_alb = qsca / qext
    
    volume_water_sphere = np.pi * 4 / 3 * ((radius)**3) # um3
    mass_water_sphere = volume_water_sphere * (1000 * 10**(-18)) # m3 to um3
    ext_cff_mss = np.pi * (radius * 10**(-6))**2 * qext / mass_water_sphere 
    sca_cff_vlm = np.pi * (radius)**2 * qsca  / volume_water_sphere * 1e6
    abs_cff_vlm = np.pi * (radius)**2 * qabs / volume_water_sphere * 1e6
    ext_cff_vlm = np.pi * (radius)**2 * qext / volume_water_sphere * 1e6

    return assym, ss_alb, sca_cff_vlm, abs_cff_vlm, ext_cff_vlm, ext_cff_mss


def net_cdf_out(asymmetry, 
                ssa, 
                sca_cff_vlm, 
                abs_cff_vlm, 
                ext_cff_vlm, 
                ext_cff_mss,
                radius,
                path_to_save, 
                medium_type,
                description,
                filename
                ):

    file = pd.DataFrame()
    file["asm_prm"] = asymmetry
    file["ss_alb"] = ssa
    file["ext_cff_mss"] = ext_cff_mss # unit m2 m-3
    file["sca_cff_vlm"] = sca_cff_vlm # unit m2 m-3
    file["abs_cff_vlm"] = abs_cff_vlm # unit m2 m-3
    file["ext_cff_vlm"] = ext_cff_vlm # unit m2 m-3
    file = file.to_xarray()
    file.attrs["medium_type"] = medium_type
    file.attrs["water_radius_um"] = radius
    file.attrs["description"] = description
    file.to_netcdf(
        str(f'{path_to_save}/{filename}_{radius}.nc'))
    

    return

# --------------------------------------------------------------------------------------
# FUNCTON CALL
# --------------------------------------------------------------------------------------
base_path = "/Users/au660413/Desktop/github/biosnicar-py/Data/OP_data"
ref_index_water = pd.read_csv(
    f'{base_path}/refractive_index_water_273K.csv')
ref_index_ice = xr.open_dataset(
    f'{base_path}/480band/rfidx_ice.nc')
path_to_save = f'{base_path}/480band/water_spherical_grains_in_air'
n_ice = ref_index_ice.re_Pic16.values
n_water = ref_index_water.n.values
k_water = ref_index_water.k.values
bbl_files = sorted(
    glob.glob(
        f'{base_path}/480band/bubbly_ice_files/*.nc'))
bbl_sizes = [int(f.split('_')[-1][:-3]) for f in bbl_files]
rad_list = list(np.array(bbl_sizes)[np.array(bbl_sizes) % 10 == 0])

n_in = n_water
k_in = k_water
n_ext = 1
medium_type = 'air'
description = "Optical properties of spherical water grains in air"
file_name = 'water_grain_in_air'

for radius in rad_list:
    assym, ssa, sca_cff_vlm, abs_cff_vlm, ext_cff_vlm, ext_cff_mss = calc_mie_params(radius, n_in, k_in, n_ext)
    net_cdf_out(assym, 
                ssa, 
                sca_cff_vlm, 
                abs_cff_vlm, 
                ext_cff_vlm,
                ext_cff_mss, 
                radius,
                path_to_save,
                medium_type, 
                description, 
                file_name
                )
    print(f'RADIUS {radius} done')

