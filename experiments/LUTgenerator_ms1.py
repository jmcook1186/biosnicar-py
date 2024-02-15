#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 18:42:05 2021

@author: Lou-Anne Chevrollier, University of Aarhus
"""

import sys 
sys.path.append("./src")
import itertools
import multiprocessing as mp
import pandas as pd
import numpy as np
import math as m
import time 


def call_SNICAR(z, density, grain_size, sza, lwc, alg): 
    
    import collections as c
    from pathlib import Path
    from snicar_feeder import snicar_feeder

    # --------------------------------------------------------------------------------------
    # 1) Initialize Inputs of the model
    # --------------------------------------------------------------------------------------
    
    Inputs = c.namedtuple(
        "Inputs",
        [
            "dir_base",
            "verbosity",
            "rf_ice",
            "incoming_i",
            "direct",
            "layer_type",
            "cdom_layer",
            "aprx_typ",
            "delta",
            "solzen",
            "toon",
            "add_double",
            "R_sfc",
            "dz",
            "rho_layers",
            "lwc",
            "lwc_pct_bubbles",
            "grain_rds",
            "side_length",
            "depth",
            "rwater",
            "nbr_lyr",
            "nbr_aer",
            "grain_shp",
            "shp_fctr",
            "grain_ar",
            "GA_units",
            "SA_units",
            "c_factor_GA",
            "c_factor_SA",
            "mss_cnc_soot1",
            "mss_cnc_soot2",
            "mss_cnc_brwnC1",
            "mss_cnc_brwnC2",
            "mss_cnc_dust1",
            "mss_cnc_dust2",
            "mss_cnc_dust3",
            "mss_cnc_dust4",
            "mss_cnc_dust5",
            "mss_cnc_ash1",
            "mss_cnc_ash2",
            "mss_cnc_ash3",
            "mss_cnc_ash4",
            "mss_cnc_ash5",
            "mss_cnc_ash_st_helens",
            "mss_cnc_Skiles_dust1",
            "mss_cnc_Skiles_dust2",
            "mss_cnc_Skiles_dust3",
            "mss_cnc_Skiles_dust4",
            "mss_cnc_Skiles_dust5",
            "mss_cnc_GreenlandCentral1",
            "mss_cnc_GreenlandCentral2",
            "mss_cnc_GreenlandCentral3",
            "mss_cnc_GreenlandCentral4",
            "mss_cnc_GreenlandCentral5",
            "mss_cnc_Cook_Greenland_dust_L",
            "mss_cnc_Cook_Greenland_dust_C",
            "mss_cnc_Cook_Greenland_dust_H",
            "mss_cnc_snw_alg",
            "mss_cnc_glacier_algae",
            "file_soot1",
            "file_soot2",
            "file_brwnC1",
            "file_brwnC2",
            "file_dust1",
            "file_dust2",
            "file_dust3",
            "file_dust4",
            "file_dust5",
            "file_ash1",
            "file_ash2",
            "file_ash3",
            "file_ash4",
            "file_ash5",
            "file_ash_st_helens",
            "file_Skiles_dust1",
            "file_Skiles_dust2",
            "file_Skiles_dust3",
            "file_Skiles_dust4",
            "file_Skiles_dust5",
            "file_GreenlandCentral1",
            "file_GreenlandCentral2",
            "file_GreenlandCentral3",
            "file_GreenlandCentral4",
            "file_GreenlandCentral5",
            "file_Cook_Greenland_dust_L",
            "file_Cook_Greenland_dust_C",
            "file_Cook_Greenland_dust_H",
            "file_snw_alg",
            "file_glacier_algae",
            "tau",
            "g",
            "SSA",
            "mu_not",
            "nbr_wvl",
            "wvl",
            "Fs",
            "Fd",
            "L_snw",
            "flx_slr",
        ],
        )


    Inputs.dir_base = str(Path(__file__).parent.parent) + "/"
    Inputs.direct = int(isinstance(sza,int))  # 1 = direct-beam, 0 = Diffuse flux
    Inputs.aprx_typ = 1  # 1 = Eddington, 2 = Quadrature, 3 = Hemispheric Mean
    Inputs.delta = 1  # 1 = Apply delta approximation, 0 = No delta
    if Inputs.direct == 1:
    	Inputs.solzen = sza  # solar zenith angle between 0 (nadir) and 89 (horizon)
    if Inputs.direct == 0:
    	Inputs.solzen = 50  # solar zenith angle between 0 (nadir) and 89 (horizon)
    # CHOOSE ATMOSPHERIC PROfile for surface-incident flux:
    #    0 = mid-latitude winter
    #    1 = mid-latitude summer
    #    2 = sub-Arctic winter
    #    3 = sub-Arctic summer
    #    4 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
    #    5 = High Mountain (summer, surface pressure of 556 hPa)
    #    6 = Top-of-atmosphere
    # NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
    # and cloudy-sky spectral fluxes are loaded when direct_beam=0
    Inputs.incoming_i = 1

    # --------------------------------------------------------------------------------------
    # 4) SET UP ICE/SNOW LAYERS
    # For granular layers only, choose toon
    # For granular layers + Fresnel layers, choose add_double
    # --------------------------------------------------------------------------------------
    
    Inputs.toon = False  # toggle toon et al tridiagonal matrix solver
    Inputs.add_double = True  # toggle addintg-doubling solver
    Inputs.dz = [0.02, z] # thickness of each vertical layer (unit = m)
    Inputs.nbr_lyr = len(Inputs.dz)  # number of snow layers
    Inputs.cdom_layer = [0]*len(Inputs.dz)  # Only for layer type == 1
    Inputs.rho_layers = [density]*len(Inputs.dz) # density of each layer (unit = kg m-3)
    Inputs.layer_type = [4, 4]  # 0 = ice grain layer, 
                                   # 1 = solid bubbly ice w/ fresnel layer,  
                                   # 2 = liquid water film, 
                                   # 3 = solid bubbly ice w/out fresnel layer
                                   # 4 = mixed water and ice grains
    Inputs.lwc = [lwc]*len(Inputs.dz)
    Inputs.lwc_pct_bubbles = 1 # amount of water as bubbles
    Inputs.nbr_wvl = 480
    Inputs.R_sfc = np.genfromtxt(
        Inputs.dir_base + "Data/OP_data/480band/r_sfc/blue_ice_spectrum_s10290721.csv",
        delimiter="csv",
    )
    Inputs.R_sfc = np.array([0 for i in range(Inputs.nbr_wvl)])


    # define source of ice refractive index data.
    # 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
    Inputs.rf_ice = 2
    
    # Ice grain shape can be
    # 0 = sphere,
    # 1 = spheroid,
    # 2 = hexagonal plate,
    # 3 = koch snowflake,
    # 4 = hexagonal prisms
    
    Inputs.grain_shp = [0]*len(Inputs.dz)  # grain shape (He et al. 2016, 2017)
    Inputs.grain_rds = [grain_size]*len(Inputs.dz) # effective grain or bubble radius
    Inputs.rwater = [0]*len(Inputs.dz) # int(Inputs.grain_rds[1]*1.055)]  # radius of optional liquid water coating

    # For 4:
    Inputs.side_length = [10000]*len(Inputs.dz)
    Inputs.depth = [10000]*len(Inputs.dz)
    
    # Shape factor = ratio of nonspherical grain effective
    # radii to that of equal-volume sphere
    # only activated when sno_shp > 1 (i.e. nonspherical)
    # 0=use recommended default value (He et al. 2017)
    # use user-specified value (between 0 and 1)
    Inputs.shp_fctr = [0]*len(Inputs.dz)
    
    # Aspect ratio (ratio of width to length)
    Inputs.grain_ar = [0]*len(Inputs.dz)
    
    # --------------------------------------------------------------------------------------
    # 5) SET LAP CHARACTERISTICS
    # --------------------------------------------------------------------------------------
    
    # Define total number of different LAPs/aerosols in model
    Inputs.nbr_aer = 30
    
    # define units for algae absorption cross section input file
    # 0 = m2/kg for MAC, ppb for mss_cnc (this is default)
    # 1 = m2/cell for MAC, cell/mL for mss_cnc
    Inputs.GA_units = 1  # glacier algae
    Inputs.SA_units = 1  # snow algae

    # determine C_factor (can be None or a number)
    # this is the concentrating factor that accounts for
    # resolution difference in field samples and model layers
    Inputs.c_factor_GA = 0 # int(0.02 / lwfilm_depth)
    Inputs.c_factor_SA = 0

    # Set names of files containing the optical properties of these LAPs:
    # uncoated BC (Bohren and Huffman, 1983)
    Inputs.file_soot1 = "mie_sot_ChC90_dns_1317.nc"
    # coated BC (Bohren and Huffman, 1983)
    Inputs.file_soot2 = "miecot_slfsot_ChC90_dns_1317.nc"
    # uncoated brown carbon (Kirchstetter et al. (2004).)
    Inputs.file_brwnC1 = "brC_Kirch_BCsd.nc"
    # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    Inputs.file_brwnC2 = "brC_Kirch_BCsd_slfcot.nc"
    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    Inputs.file_dust1 = "dust_balkanski_central_size1.nc"
    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    Inputs.file_dust2 = "dust_balkanski_central_size2.nc"
    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    Inputs.file_dust3 = "dust_balkanski_central_size3.nc"
    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    Inputs.file_dust4 = "dust_balkanski_central_size4.nc"
    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    Inputs.file_dust5 = "dust_balkanski_central_size5.nc"
    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    Inputs.file_ash1 = "volc_ash_eyja_central_size1.nc"
    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    Inputs.file_ash2 = "volc_ash_eyja_central_size2.nc"
    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    Inputs.file_ash3 = "volc_ash_eyja_central_size3.nc"
    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    Inputs.file_ash4 = "volc_ash_eyja_central_size4.nc"
    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    Inputs.file_ash5 = "volc_ash_eyja_central_size5.nc"
    # ashes from Mount Saint Helen's
    Inputs.file_ash_st_helens = "volc_ash_mtsthelens_20081011.nc"
    # Colorado dust 1 (Skiles et al 2017)
    Inputs.file_Skiles_dust1 = "dust_skiles_size1.nc"
    # Colorado dust 2 (Skiles et al 2017)
    Inputs.file_Skiles_dust2 = "dust_skiles_size2.nc"
    # Colorado dust 3 (Skiles et al 2017)
    Inputs.file_Skiles_dust3 = "dust_skiles_size3.nc"
    # Colorado dust 4 (Skiles et al 2017)
    Inputs.file_Skiles_dust4 = "dust_skiles_size4.nc"
    # Colorado dust 5 (Skiles et al 2017)
    Inputs.file_Skiles_dust5 = "dust_skiles_size5.nc"
    # Greenland dust 1 (Polashenski et al 2015)
    Inputs.file_GreenlandCentral1 = "dust_greenland_central_size1.nc"
    # Greenland dust 2 (Polashenski et al 2015)
    Inputs.file_GreenlandCentral2 = "dust_greenland_central_size2.nc"
    # Greenland dust 3 (Polashenski et al 2015)
    Inputs.file_GreenlandCentral3 = "dust_greenland_central_size3.nc"
    # Greenland dust 4 (Polashenski et al 2015)
    Inputs.file_GreenlandCentral4 = "dust_greenland_central_size4.nc"
    # Greenlanddust 5 (Polashenski et al 2015)
    Inputs.file_GreenlandCentral5 = "dust_greenland_central_size5.nc"
    # Dark Zone dust 1 (Cook et al. 2019 "LOW") NOT FUNCTIONAL (COMING SOON)
    Inputs.file_Cook_Greenland_dust_L = "dust_greenland_Cook_LOW_20190911.nc"
    # Dark Zone dust 2 (Cook et al. 2019 "mean") NOT FUNCTIONAL (COMING SOON)
    Inputs.file_Cook_Greenland_dust_C = "dust_greenland_Cook_CENTRAL_20190911.nc"
    # Dark Zone dust 3 (Cook et al. 2019 "HIGH") NOT FUNCTIONAL (COMING SOON)
    Inputs.file_Cook_Greenland_dust_H = "dust_greenland_Cook_HIGH_20190911.nc"
    # Snow Algae (Cook et al. 2017a, spherical, C nivalis)
    Inputs.file_snw_alg = "snow_algae_empirical_Chevrollier2023.nc"
    # Glacier Algae (Cook et al. 2020)
    Inputs.file_glacier_algae =  "GA_Chevrollier2022_r4.9_L18.8.nc"

    # Indicate mass mixing ratios scenarios for each impurity
    # default units are ng(species)/g(ice), or ppb
    # However, snow and glacier algae can be provided in in cells/mL.
    # To use cells/mL for algae, set GA_units == 1.
    # The script will loop over the different mixing scenarios
    
    Inputs.mss_cnc_soot1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_soot2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_brwnC1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_brwnC2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_dust1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_dust2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_dust3 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_dust4 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_dust5 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash3 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash4 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash5 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_ash_st_helens = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Skiles_dust1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Skiles_dust2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Skiles_dust3 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Skiles_dust4 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Skiles_dust5 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_GreenlandCentral1 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_GreenlandCentral2 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_GreenlandCentral3 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_GreenlandCentral4 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_GreenlandCentral5 = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Cook_Greenland_dust_L = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Cook_Greenland_dust_C = [0] * len(Inputs.dz)
    Inputs.mss_cnc_Cook_Greenland_dust_H = [0] * len(Inputs.dz)
    Inputs.mss_cnc_glacier_algae = [0] * len(Inputs.dz)
    Inputs.mss_cnc_snw_alg = [alg, 0] 



    # --------------------------------------------------------------------------------------
    # CALL FUNCTIONS AND PLOT OUTPUTS
    # --------------------------------------------------------------------------------------


    Outputs = snicar_feeder(Inputs)
    albedo = Outputs.albedo
    BBA = Outputs.BBA
    # data = np.float16(np.append(albedo[70:170], BBA))
    data = np.float16(np.append(albedo, BBA))

    return data


#############################################################################
################ INPUTS
#############################################################################
path_to_save_files = '/data/lou'

## PARAMS FOR CHLA ABS FEATURE 
densities = [500+i*50 for i in range(0,5)] # maybe this can even be removed
lwcs = [0.05] # maybe can be fixed or rm
sza_list = [40, 45, 50, 55, 60] # maybe can be removed, use diff instead
algs = [1e3, 5e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 
        1e5, 1.2e5,  1.4e5, 1.6e5, 1.8e5,  
        2e5, 2.25e5, 2.5e5, 2.75e5, 3e5, 3.25e5, 3.5e5, 3.75e5, 4e5, 4.25e5,
        4.5e5, 4.75e5, 5e5, 5.5e5, 6e5, 6.5e5, 7e5, 7.5e5, 8e5, 8.5e5, 9e5, 9.5e5,
        1e6, 1.1e6, 1.2e6, 1.3e6, 1.4e6, 1.5e6, 1.6e6, 1.7e6, 1.8e6, 1.9e6,
        2e6, 2.2e6, 2.4e6, 2.6e6, 2.8e6, 3e6, 3.25e6, 3.5e6, 3.75e6, 4e6,
        4.5e6, 5e6, 5.5e6, 6e6, 6.5e6, 7e6, 7.5e6, 8e6, 9e6, 1e7,
        1.1e7, 1.2e7, 1.3e7, 1.4e7, 1.6e7, 1.8e7, 2.1e7, 2.5e7, 3.1e7, 4e7,
        5e7, 6.5e7, 1e8, 2e8, 1e9]
grain_sizes = [500 + i * 40 for i in range(64)] 
depths = [1]

## PARAMS FOR RED SPECTRA NIR INVERSION
# densities = [500+i*20 for i in range(0,11)]
# lwcs = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2] 
# depths = [1]
# sza_list = [43, 42, 44, 40, 48, 'diff']  # 39, 46, 47
# grain_sizes = [500, 550, 600,
#              650, 700, 750,
#              800, 900, 1000, 1100,
#              1200, 1300, 1400, 1500,
#              2000, 3000]
# algs = [0,2500, 5000, 7500, 10000,
#         12500,
#         15000, 20000, 25000,
#         30000,35000, 
#         40000,45000, 
#         50000,55000,
#         60000,65000, 
#         70000,75000,
#         80000,85000,
#         90000,95000,
#         100000]


# #### LUT PARAMETERS FIELD SPECTRA
# densities = [
# 380, 
# 400, 
# 420, 
# 440, 
# 460, 
# 480, 
# 500,
# 520, 
# 540, 
# 560, 
# 580, 
# 600, 
# 620, 
# 640, 
# 660, 
# 680, 
# 700, 
# 720, 
# 740, 
# 760, 
# 780, 
# 800, 
# 820, 
# 840, 
# 860, 
# 880,
# 900, 
# 916.9
# ] # 28

# sza_list = [42, 
# 'diff', 
# 52, 
# 54, 
# 46, 
# 43, 
# 44, 
# 47, 
# 45] # 9

# depths = [
# 5e-5, 
# 3e-4,
# 6.3e-4,
# 1e-3,
# 1.5e-3,
# 2e-3,
# 2.5e-3,
# 3e-3,
# 3.5e-3,
# 4.2e-3,
# 5e-3, 
# 6e-3,
# 8e-3,
# 1.1e-2,
# 1] # 15

# lwcs = [0,
# 0.02,
# 0.04,
# 0.06,
# 0.08,
# 0.1,
# 0.12,
# 0.14,
# 0.16,
# 0.18,
# 0.2,
# 0.22,
# 0.24,
# 0.26,
# 0.28,
# 0.3] # 16

# grain_sizes = [
# [50],
# [55],
# [60],
# [65],
# [70],
# [75],
# [80],
# [85],
# [90],
# [95],
# [100],
# [110],
# [120],
# [130],
# [140],
# [150],
# [160],
# [170],
# [180],
# [190],
# [200],
# [210],
# [220],
# [230],
# [240],
# [250],
# [260],
# [270],
# [280],
# [290],
# [300],
# [325],
# [350],
# [375],
# [400],
# [425],
# [450],
# [475],
# [500],
# [550],
# [600],
# [650],
# [700],
# [750],
# [800],
# [850],
# [900],
# [950],
# [1000],
# [1100],
# [1200],
# [1300],
# [1400],
# [1500],
# [1600],
# [1700],
# [1800],
# [1900],
# [2000],
# [2100],
# [2200],
# [2300],
# [2500],
# [2750],
# [3000],
# [3250],
# [3500],
# [4000],
# [4500],
# [5000],
# [6000],
# [7000],
# [8000],
# [9000],
# [11000],
# [13000],
# [16000],
# [20000]
# ] # 78

# #### LUT PARAMETERS PRISMA
# densities = [
# 380, 
# 400, 
# 420, 
# 440, 
# 460, 
# 480, 
# 500,
# 520, 
# 540, 
# 560, 
# 580, 
# 600, 
# 620, 
# 640, 
# 660, 
# 680, 
# 700, 
# 720, 
# 740, 
# 760, 
# 780, 
# 800, 
# 820, 
# 840, 
# 860, 
# 880,
# 900, 
# 916.9
# ] # 28

# sza_list = [47]

# depths = [
# 5e-5, 
# 3e-4,
# 6.3e-4,
# 1e-3,
# 1.5e-3,
# 2e-3,
# 2.5e-3,
# 3e-3,
# 3.5e-3,
# 4.2e-3,
# 5e-3, 
# 6e-3,
# 8e-3,
# 1.1e-2,
# 1] # 15

# lwcs = [
# [0],
# [0.02],
# [0.04],
# [0.06],
# [0.08],
# [0.1],
# [0.12],
# [0.14],
# [0.16],
# [0.18],
# [0.2],
# [0.22],
# [0.24],
# [0.26],
# [0.28],
# [0.3]
# ] # 16

# grain_sizes = [
# 50,
# 55,
# 60,
# 65,
# 70,
# 75,
# 80,
# 85,
# 90,
# 95,
# 100,
# 110,
# 120,
# 130,
# 140,
# 150,
# 160,
# 170,
# 180,
# 190,
# 200,
# 210,
# 220,
# 230,
# 240,
# 250,
# 260,
# 270,
# 280,
# 290,
# 300,
# 310,
# 320,
# 330,
# 340,
# 350,
# 360,
# 370,
# 380,
# 390,
# 400,
# 410,
# 420,
# 430,
# 440,
# 450,
# 460,
# 470,
# 480,
# 490,
# 500,
# 520,
# 540,
# 560,
# 580,
# 600,
# 620,
# 640,
# 660,
# 680,
# 700,
# 720,
# 740,
# 760,
# 780,
# 800,
# 820,
# 840,
# 860,
# 880,
# 900,
# 920,
# 940,
# 960,
# 980,
# 1000,
# 1020,
# 1040,
# 1060,
# 1080,
# 1100,
# 1120,
# 1140,
# 1160,
# 1180,
# 1200,
# 1220,
# 1240,
# 1260,
# 1280,
# 1300,
# 1320,
# 1340,
# 1360,
# 1380,
# 1400,
# 1420,
# 1440,
# 1460,
# 1480,
# 1500,
# 1550,
# 1600,
# 1650,
# 1700,
# 1750,
# 1800,
# 1850,
# 1900,
# 1950,
# 2000,
# 2050,
# 2100,
# 2150,
# 2200,
# 2250,
# 2300,
# 2350,
# 2400,
# 2450,
# 2500,
# 2550,
# 2600,
# 2650,
# 2700,
# 2750,
# 2800,
# 2850,
# 2900,
# 2950,
# 3000,    
# 3250,
# 3500,  
# 4000,
# 4500,
# 5000,
# 5500,
# 6000,
# 6500,
# 7000,
# 7500,
# 8000,
# 8500,
# 9000,
# 10000,
# 20000
# ]



#############################################################################
################ CALL MODEL
#############################################################################
start = time.time()
paramlist = list(set((itertools.product(depths, densities, grain_sizes, 
                                          sza_list, lwcs, algs))))
if __name__ == '__main__':
    nb_cores = 150
    pool = mp.Pool(nb_cores)
    print(f'starting simulation on {nb_cores} cores')
    data = pool.starmap(call_SNICAR, paramlist)
    # pool.close()  
    # pool.join() 
    df = pd.DataFrame.from_records(data)
    df = df.transpose()
    df.columns = [str(i) for i in paramlist]
    df.to_feather(f'{path_to_save_files}/021524_lut_snow_with_snw_alg.feather', 
                  compression='zstd')
    print('time for 1 lut: {}'.format(time.time() - start))


# for lwc in lwcs:
#  	start = time.time()
#  	paramlist = list(set((itertools.product(depths,  
#                                             densities, 
#                                             grain_sizes, 
#                                             sza_list, 
#                                             lwc))))
#  	if __name__ == '__main__':
#           nb_cores = 200
#           pool = mp.Pool(nb_cores)
#           print(f'starting simulation on {nb_cores} cores')
#           data = pool.starmap(call_SNICAR,paramlist)
#           pool.close()  
#           pool.join()
#           df = pd.DataFrame.from_records(data)
#           df = df.transpose()
#           df.columns = [str(i) for i in paramlist]
#           df.to_feather(f'{path_to_save_files}/020523_lut_{lwc[0]}_lwc_for_prisma.feather', 
#                 compression='zstd') 
#           print('time for 1 lut: {}'.format(time.time() - start))
         
# start = time.time()
# paramlist = list(set((itertools.product(depths, lwfilm_dz, densities, grain_sizes, sza_list, algae))))

# if __name__ == '__main__':
#     nb_cores = 80
#     pool = mp.Pool(nb_cores)
#     print(f'starting simulation on {nb_cores} cores')
#     data = pool.starmap(call_SNICAR,paramlist)
#     pool.close()  
#     pool.join() 
#     df = pd.DataFrame.from_records(data)
#     df = df.transpose()
#     df.columns = [str(i) for i in paramlist]
#     df.to_feather(f'{path_to_save_files}/271023_lut_snow_with_sa.feather', 
#                   compression='zstd')
#     print('time for 1 lut: {}'.format(time.time() - start))

## 20.15mn on 100 cores for 136080 runs casted to float 16 -- size 128MB
## 17.8mn on 100 cores for 136080 runs not casted -- size 228MB
## 21 mn on 100 cores for 136080 runs to float 32 -- size 162 MB
## 27 mn on 100 cores for 136080 runs casted to float 16 but full size -- size 248MB


