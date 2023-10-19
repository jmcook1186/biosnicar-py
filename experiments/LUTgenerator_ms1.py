#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 18:42:05 2021

@author: au660413
"""
import sys 
sys.path.append("./src")
import itertools
import multiprocessing as mp
import pandas as pd
import numpy as np
import math as m
import time 


def call_SNICAR(z, lwfilm_dz, density, grain_size, sza, lwc): 
    
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
    Inputs.dz = [lwfilm_dz, z, 1] # thickness of each vertical layer (unit = m)
    Inputs.nbr_lyr = len(Inputs.dz)  # number of snow layers
    Inputs.cdom_layer = [0]*len(Inputs.dz)  # Only for layer type == 1
    Inputs.rho_layers = [916.999,density, density]  # density of each layer (unit = kg m-3)
    Inputs.layer_type = [2, 3, 3]  # 0 = ice grain layer, 
                                   # 1 = solid bubbly ice w/ fresnel layer,  
                                   # 2 = liquid water film, 
                                   # 3 = solid bubbly ice w/out fresnel layer
    Inputs.lwc = lwc
    Inputs.lwc_type = 1 # 0 = optically equivalent lwc, 
                        # 1 = water bubbles w/ same size as air bubbles
    Inputs.nbr_wvl = 480
    Inputs.R_sfc = np.genfromtxt(
        Inputs.dir_base + "Data/OP_data/480band/r_sfc/blue_ice_spectrum_s10290721.csv",
        delimiter="csv",
    )
    # Inputs.R_sfc = np.array([0 for i in range(Inputs.nbr_wvl)])


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
    Inputs.grain_rds = [grain_size, grain_size, grain_size] #[20000, grain_size]  # effective grain or bubble radius
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
    Inputs.file_snw_alg = "SA_Chevrollier2022_r8.99.nc"
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
    Inputs.mss_cnc_snw_alg = [0] * len(Inputs.dz)
    Inputs.mss_cnc_glacier_algae = [0] * len(Inputs.dz)


    # --------------------------------------------------------------------------------------
    # CALL FUNCTIONS AND PLOT OUTPUTS
    # --------------------------------------------------------------------------------------


    Outputs = snicar_feeder(Inputs)
    # if Inputs.layer_type[1] == 1:
        # SSA = (6* (917 - Inputs.rho_layers[1]) / 917 / 
       # (Inputs.rho_layers[1] * Inputs.grain_rds[1] * 2 * 10 ** (-6)))
    # if Inputs.layer_type[1] == 0:
        # SSA = 3/(Inputs.rho_layers[1] * (Inputs.grain_rds[1]* 10 ** (-6)))
    albedo = Outputs.albedo
    BBA = Outputs.BBA
    data = np.append(albedo, BBA)
    
    #calculate band ratios using spectral responses from S3
    #signal calculated at each wvl by linear interp and 
    #multiplied by weight and finally summed 
    #srf = pd.read_csv('experiments/srf_resampled.csv')
    #albedo_interp = pd.DataFrame().reindex_like(srf)
    #bands = np.zeros(21)
    #for band in range(1,22): 
    #     wvl_band = srf[f'wvl{band}']
    #     albedo_interp[f'wvl{band}']=srf[f'wvl{band}']
         # systematically find two values of albedo to interpolate
    #     for i,val in enumerate(wvl_band):
    #           if not m.isnan(val):
    #           #if the wavelength finishes by 5, can be picked from albedo directly 
    #                if (val%5==0 and val%10!=0):
    #                    albedo_interp.iloc[i][f'band{band}'] = albedo[int((val-5)/10-20)]
    #           #if the wavelength doesnt, find closest values to interp:
    #                
    #                else:
    #                    wvl_up=0
    #                    wvl_down=0
    #                    if  val%10==0:
    #                        wvl_down = val - 5
    #                        wvl_up = val + 5
    #                    elif (val - val%5)%10 == 0:
    #                        wvl_up = val + 5 - val%5
    #                        wvl_down = wvl_up - 10
    #                    else:
    #                        wvl_down = val - val%5
    #                        wvl_up = wvl_down + 10
    #                    new_alb = albedo[int((wvl_down-5)/10-20)] +  (val-wvl_down) * (albedo[int((wvl_up-5)/10-20)]-albedo[int((wvl_down-5)/10-20)]) / (wvl_up-wvl_down)
    #                    albedo_interp.iloc[i][f'band{band}'] = new_alb
    #     bands[band-1]=np.nansum(np.array(albedo_interp[f'band{band}'])*np.array(srf[f'band{band}']))              
    #data = np.append(data,bands)
    #data = np.append(data,ga)
    return data.astype(np.float16)



# densities=[300+i*40 for i in range(0,16)] # [300+i*20 for i in range(0,21)]
# grain_sizes= [1000, 2000, 3000,4000,5000,6000,7000,8000,
# 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000]
# #grain_sizes=[1500,2000, 2500,3000,
# #3500, 4000, 4500, 5000,5500,6000,6500,7000,7500,8000,
# #8500, 9000,9500, 10000, 10500, 11000, 11500, 12000, 12500,13000, 13500,14000,
# #14500, 15000,15500,16000,16500,17000, 17500, 18000, 18500, 19000, 19500, 20000]

# dz = [round(elem, 3) 
# for elem in np.concatenate((np.arange(0.01,0.1, 0.002), np.arange(0.1,0.15, 0.005), np.arange(0.15, 0.31, 0.01)))]

# dz1 = [round(elem, 3) 
# 	for elem in np.arange(0.03,0.6, 0.02)]
# # dz1 = [round(elem, 3) for elem in np.concatenate((np.arange(0.02,0.1, 0.02), np.arange(0.1, 1.05, 0.1)))]

# #dz2 = [round(elem, 3) 
# #	for elem in np.concatenate((np.arange(0.001,0.005, 0.001),np.arange(0.01,1, 0.025)))]
# # remove values that produce SSA too high
# SSA = ({(density, grain_size): (6* (917 - density) / 917 / 
#        (density * grain_size * 2 * 10 ** (-6))) 
#         for density in densities for grain_size in grain_sizes})
# SSA = ({(density, grain_size): 3/(density * (grain_size* 10 ** (-6)))
# 	for density in densities for grain_size in grain_sizes})
# params_to_remove = [] #k for k in SSA.keys() if SSA[k] > 4]
# algae = [0,5000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,120000,140000,150000]
# sza_list = [[42], [44], [46], [48], [50], [52], [54], [56], [58], [60]]

#############################################################################
################ INPUTS
#############################################################################
path_to_save_files = '/data/lou'
## PARAMS FOR EQUIVALENT/WATER SPHERES IN ICE LUTS
densities = [330+i*10 for i in range(0,59)]
grain_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900,
                1000, 2000, 3000,
                4000, 5000,  6000, 
                7000, 8000, 9000, 
                10000, 11000, 12000,  
                13000, 14000, 15000, 
                16000, 17000, 18000, 
                19000, 20000] 
sza_list = [[43], [42], 
            [46], [48], 
            [52], [44], 
            [47], [54],
            [49], ['diff']] 
depths = [5e-5, 1e-4,  
          2e-4, 2.5e-4, 
          3e-4, 4e-4, 
          5e-4, 1e-3,
          2.5e-3, 
          4e-3, 5e-3, 
          6e-3,7e-3,
          8e-3, 
          1e-2, 1.5e-2, 
          3e-2, 4.5e-2, 
          6e-2, 0.1,
          1] 
lwfilm_dz = [1e-10] #[1e-10, 0.0001, 0.0005]
lwcs = [[0], [0.01], [0.02], [0.03], [0.04], [0.05], [0.06], 
        [0.07], [0.08], [0.09], [0.1], [0.11], [0.12], [0.13], [0.14],  
        [0.15],[0.175], [0.2], [0.225], 
        [0.25], [0.275], [0.23]] 

## PARAMS FOR MIXED SPHERES ICE WATER 
# densities = [300+i*10 for i in range(0,77)]
# grain_sizes = [100, 200, 300, 400, 500, 600, 700, 800, 900,
#                 1000, 2000, 3000,
#                 4000, 5000,  6000, 
#                 7000, 8000, 9000, 
#                 10000, 11000, 12000,  
#                 13000, 14000, 15000, 
#                 16000, 17000, 18000, 
#                 19000, 20000] 
# depths = [1]
# lwfilm_dz = [1e-10]
# lwcs = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 
#         0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,  
#         0.15,0.175, 0.2, 0.225, 
#         0.25, 0.275, 0.23] 
# sza_list = [[43], [46], 
#             ['diff'], [48], 
#             [49], [44], 
#             [47], [52],
#             [54], [42]] 



#############################################################################
################ CALL MODEL
#############################################################################

for lwc in lwcs:
 	start = time.time()
 	paramlist = sorted(set((itertools.product(depths, lwfilm_dz, densities, grain_sizes, sza_list, lwc))))
 	if __name__ == '__main__':
         nb_cores = 80
         pool = mp.Pool(nb_cores)
         print(f'starting simulation on {nb_cores} cores')
         data = pool.starmap(call_SNICAR,paramlist)
         pool.close()  
         pool.join()
         df = pd.DataFrame.from_records(data)
         df = df.transpose()
         df.columns = [str(i) for i in paramlist]
         df.to_feather(f'{path_to_save_files}/191023_lut_{lwc}_lwc_water_bbl_in_ice_varying_fresnel_depth.feather', 
                compression='zstd') 
         print('time for 1 lut: {}'.format(time.time() - start))
# for sza in sza_list:
# 	start = time.time()
# 	paramlist = sorted(set((itertools.product(depths, lwfilm_dz, densities, grain_sizes, sza, lwcs))))
# 	if __name__ == '__main__':
# 		nb_cores = 100
# 		pool = mp.Pool(nb_cores)
# 		print(f'starting simulation on {nb_cores} cores')
# 		data = pool.starmap(call_SNICAR,paramlist)
# 		pool.close()  
# 		pool.join() 
# 		df = pd.DataFrame.from_records(data)
# 		df = df.transpose()
# 		df.columns = [str(i) for i in paramlist]
# 		df.to_feather(f'{path_to_save_files}/131023_lut_{sza[0]}_mixed_water_ice_spheres.feather', 
#                 compression='zstd')
# 		print('time for 1 lut: {}'.format(time.time() - start))


