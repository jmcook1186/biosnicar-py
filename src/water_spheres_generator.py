#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 22:38:08 2023

@author: Lou-Anne Chevrollier, University of Aarhus
"""

import pandas as pd
import numpy as np
from miepython import mie, ez_mie
import xarray as xr
import glob
import PyMieScatt

def calc_mie_params(radius, n_in, k_in, n_ext):
    wvl = np.arange(0.205, 5, 0.01) # um
    # X = 2 * np.pi * radius / wvl  # unitless
    qext, qsca, qback, g = ez_mie(n_in - 1j * k_in, 
                                  2 * radius,
                                  wvl,
                                  n_ext) 
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

def calc_mie_coated_params(radius_in, radius_total, 
                           n_in, k_in, 
                           n_shell, k_shell,
                           n_ext, wvl):
    mCore = n_in + 1j * k_in
    mShell = n_shell + 1j * k_shell
    dCore = radius_in * 1000 # um to nm
    dShell = radius_total * 1000
    nMedium = n_ext
    qext, qsca, qabs, g, qpr, qback, qratio = PyMieScatt.MieQCoreShell(
        mCore, mShell, wvl, dCore, dShell, 
        nMedium, asDict=False, asCrossSection=False)
    
    assym = g
    ss_alb = qsca / qext
    
    volume_water_sphere = np.pi * 4 / 3 * ((radius_total)**3) # um3
    mass_water_sphere = volume_water_sphere * (1000 * 10**(-18)) # m3 to um3
    ext_cff_mss = np.pi * (radius_total * 10**(-6))**2 * qext / mass_water_sphere 
    sca_cff_vlm = np.pi * (radius_total)**2 * qsca  / volume_water_sphere * 1e6
    abs_cff_vlm = np.pi * (radius_total)**2 * qabs / volume_water_sphere * 1e6
    ext_cff_vlm = np.pi * (radius_total)**2 * qext / volume_water_sphere * 1e6
    
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


#%%
# --------------------------------------------------------------------------------------
# FUNCTON CALL
# --------------------------------------------------------------------------------------
base_path = "/Users/au660413/Desktop/github/biosnicar-py/Data/OP_data"
ref_index_water = pd.read_csv(
    f'{base_path}/refractive_index_water_273K.csv')
ref_index_ice = xr.open_dataset(
    f'{base_path}/480band/rfidx_ice.nc')
n_ice = ref_index_ice.re_Pic16.values
k_ice = ref_index_ice.im_Pic16.values

n_water = ref_index_water.n.values
k_water = ref_index_water.k.values
bbl_files = sorted(
    glob.glob(
        f'{base_path}/480band/bubbly_ice_files/*.nc'))
bbl_sizes = [int(f.split('_')[-1][:-3]) for f in bbl_files]
rad_list = list(np.array(bbl_sizes)[np.array(bbl_sizes) % 10 == 5])

path_to_save = f'{base_path}/480band/water_spherical_grains_in_ice'
n_in = n_water # 
k_in =  k_water 
n_ext = n_ice
medium_type = 'ice'
description = "Optical properties of spherical water spheres in ice"
file_name = 'water_grain_in_ice'

for radius in rad_list:
    assym, ssa, sca_cff_vlm, abs_cff_vlm, ext_cff_vlm, ext_cff_mss = calc_mie_params(radius, n_in, k_in, n_ext)
    data_water_spheres_ = net_cdf_out(assym, 
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
    
    
#%% coated spheres

# path_to_save = f'{base_path}/480band/air_spherical_bubbles_coated_in_ice'

# n_in = n_water
# k_in = k_water
# n_ext = n_ice
# medium_type = 'ice'
# description = "Optical properties of spherical air spheres with water coating in ice"
# file_name = 'air_bubble_coated_in_ice'

# assym_total = []
# ssa_total = []
# ext_cff_mss_total = []
# sca_cff_vlm_total = []
# abs_cff_vlm_total = []
# ext_cff_vlm_total = []
# wavelengths = np.arange(205, 5000, 10)

# radius = 600
# for i in range(len(wavelengths)):
#     assym, ssa, sca_cff_vlm, abs_cff_vlm, ext_cff_vlm, ext_cff_mss = calc_mie_coated_params(radius*0.9, radius, 
#                                                                                             n_ice[i], k_ice[i], 
#                                                                                             n_water[i], k_water[i], 
#                                                                                             1.0, 
#                                                                                             wavelengths[i])
#     assym_total.append(assym)
#     ssa_total.append(ssa)
#     sca_cff_vlm_total.append(sca_cff_vlm)
#     abs_cff_vlm_total.append(abs_cff_vlm)
#     ext_cff_vlm_total.append(ext_cff_vlm)
#     ext_cff_mss_total.append(ext_cff_mss)
#     print(i)

# data_bbl = net_cdf_out(assym_total, 
#             ssa_total, 
#             sca_cff_vlm_total, 
#             abs_cff_vlm_total, 
#             ext_cff_vlm_total,
#             ext_cff_mss_total, 
#             radius,
#             path_to_save,
#             medium_type, 
#             description, 
#             file_name
#             )
# print(f'RADIUS {radius} done')

