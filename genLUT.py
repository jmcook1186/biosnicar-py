#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:24:16 2021

@author: au660413
"""

import pandas as pd
    
def call_SNICAR(z, input_concentration, density, layer_types, grain_size, sza):

    from SNICAR_feeder import snicar_feeder
    import matplotlib.pyplot as plt
    import numpy as np
    import collections as c

######################################
## 1) Initialize inputs of the model
######################################

    inputs = c.namedtuple('inputs',['dir_base', 'verbosity', \
        'rf_ice', 'incoming_i', 'DIRECT', 'layer_type', 'cdom_layer',\
        'APRX_TYP', 'DELTA', 'solzen', 'TOON', 'ADD_DOUBLE', 'R_sfc', 'dz', 'rho_layers', 'grain_rds',\
        'side_length', 'depth', 'rwater', 'nbr_lyr', 'nbr_aer', 'grain_shp', 'shp_fctr', 'grain_ar', 'GA_units', 'SA_units',\
        'Cfactor_GA','Cfactor_SA','mss_cnc_soot1', 'mss_cnc_soot2', 'mss_cnc_brwnC1', 'mss_cnc_brwnC2', 'mss_cnc_dust1',\
        'mss_cnc_dust2', 'mss_cnc_dust3', 'mss_cnc_dust4', 'mss_cnc_dust5', 'mss_cnc_ash1', 'mss_cnc_ash2',\
        'mss_cnc_ash3', 'mss_cnc_ash4', 'mss_cnc_ash5', 'mss_cnc_ash_st_helens', 'mss_cnc_Skiles_dust1', 'mss_cnc_Skiles_dust2',\
        'mss_cnc_Skiles_dust3', 'mss_cnc_Skiles_dust4', 'mss_cnc_Skiles_dust5', 'mss_cnc_GreenlandCentral1',\
        'mss_cnc_GreenlandCentral2', 'mss_cnc_GreenlandCentral3', 'mss_cnc_GreenlandCentral4',\
        'mss_cnc_GreenlandCentral5', 'mss_cnc_Cook_Greenland_dust_L', 'mss_cnc_Cook_Greenland_dust_C',\
        'mss_cnc_Cook_Greenland_dust_H', 'mss_cnc_snw_alg', 'mss_cnc_glacier_algae', 'FILE_soot1',\
        'FILE_soot2', 'FILE_brwnC1', 'FILE_brwnC2', 'FILE_dust1', 'FILE_dust2', 'FILE_dust3', 'FILE_dust4', 'FILE_dust5',\
        'FILE_ash1', 'FILE_ash2', 'FILE_ash3', 'FILE_ash4', 'FILE_ash5', 'FILE_ash_st_helens', 'FILE_Skiles_dust1', 'FILE_Skiles_dust2',\
        'FILE_Skiles_dust3', 'FILE_Skiles_dust4', 'FILE_Skiles_dust5', 'FILE_GreenlandCentral1',\
        'FILE_GreenlandCentral2', 'FILE_GreenlandCentral3', 'FILE_GreenlandCentral4', 'FILE_GreenlandCentral5',\
        'FILE_Cook_Greenland_dust_L', 'FILE_Cook_Greenland_dust_C', 'FILE_Cook_Greenland_dust_H', 'FILE_snw_alg', 'FILE_glacier_algae',\
        'tau', 'g', 'SSA', 'mu_not', 'nbr_wvl', 'wvl', 'Fs', 'Fd', 'L_snw', 'flx_slr'])
    
    
    ##############################
    ## 2) Set working directory 
    ##############################
    
    # set dir_base to the location of the BioSNICAR_GO_PY folder
    inputs.dir_base = '/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/'
    savepath = inputs.dir_base # base path for saving figures
    write_config_to_textfile = False # toggle to save config to file
    inputs.verbosity = 0 # 1 to print real-time updates
    
    
    ################################
    ## 3) Choose plot/print options
    ################################
    
    show_figs = True # toggle to display spectral albedo figure
    save_figs = False # toggle to save spectral albedo figure to file
    print_BBA = True # toggle to print broadband albedo to terminal
    print_band_ratios = False # toggle to print various band ratios to terminal
    smooth = False # apply optional smoothing function (Savitzky-Golay filter)
    window_size = 9 # if applying smoothing filter, define window size
    poly_order = 3 # if applying smoothing filter, define order of polynomial
    
    #######################################
    ## 4) RADIATIVE TRANSFER CONFIGURATION
    #######################################
    
    inputs.DIRECT   = 1       # 1= Direct-beam incident flux, 0= Diffuse incident flux
    inputs.APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    inputs.DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
    inputs.solzen   = sza      # if DIRECT give solar zenith angle between 0 and 89 degrees (from 0 = nadir, 90 = horizon)
    
    # CHOOSE ATMOSPHERIC PROFILE for surface-incident flux:
    #    0 = mid-latitude winter
    #    1 = mid-latitude summer
    #    2 = sub-Arctic winter
    #    3 = sub-Arctic summer
    #    4 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
    #    5 = High Mountain (summer, surface pressure of 556 hPa)
    #    6 = Top-of-atmosphere
    # NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
    # and cloudy-sky spectral fluxes are loaded when direct_beam=0
    inputs.incoming_i = 4
    
    ###############################################################
    ## 4) SET UP ICE/SNOW LAYERS
    # For granular layers only, choose TOON
    # For granular layers + Fresnel layers, choose ADD_DOUBLE
    ###############################################################
    
    inputs.TOON = False # toggle Toon et al tridiagonal matrix solver
    inputs.ADD_DOUBLE = True # toggle adding-doubling solver
    
    inputs.dz = z # thickness of each vertical layer (unit = m)
    inputs.nbr_lyr = len(inputs.dz)  # number of snow layers
    inputs.layer_type = layer_types # Fresnel layers for the ADD_DOUBLE option, set all to 0 for the TOON option
    inputs.cdom_layer = [0,0] # Only functional if layer type == 1, based on CDOM data from Halbach et al. 2021 (in prep)
    inputs.rho_layers = [density]*len(inputs.dz)  # density of each layer (unit = kg m-3) 
    inputs.nbr_wvl=480 
    #inputs.R_sfc = np.array([0.1 for i in range(inputs.nbr_wvl)]) # reflectance of undrlying surface - set across all wavelengths
    inputs.R_sfc = np.genfromtxt(inputs.dir_base+'/Data/rain_polished_ice_spectrum.csv', delimiter = 'csv') # import underlying ice from file
    
    ###############################################################################
    ## 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
    # For hexagonal plates or columns of any size choose GeometricOptics
    # For sphere, spheroids, koch snowflake with optional water coating choose Mie
    ###############################################################################
    
    inputs.rf_ice = 2 # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
    
    # Ice grain shape can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3 = koch snowflake, 4 = hexagonal prisms
    # For 0,1,2,3:
    inputs.grain_shp =[0,0] # grain shape(He et al. 2016, 2017)
    inputs.grain_rds = [grain_size]*len(inputs.dz)  # effective grain radius of snow/bubbly ice (becomes bubble rds when layer_type==1)
    inputs.rwater = [0, 0] # radius of optional liquid water coating
    
    # For 4:
    inputs.side_length = [800,800] 
    inputs.depth = [10000,10000]
    
    # Shape factor = ratio of nonspherical grain effective radii to that of equal-volume sphere
    ### only activated when sno_shp > 1 (i.e. nonspherical)
    ### 0=use recommended default value (He et al. 2017)
    ### use user-specified value (between 0 and 1)
    inputs.shp_fctr = [0,0] 
    
    # Aspect ratio (ratio of width to length)
    inputs.grain_ar = [0,0] 
    
    #######################################
    ## 5) SET LAP CHARACTERISTICS
    #######################################
    
    # Define total number of different LAPs/aerosols in model
    inputs.nbr_aer = 30
    
    # define units for algae absorption cross section input file
    # 0 = m2/kg for MAC, ppb for mss_cnc (this is default)
    # 1 = m2/cell for MAC, cell/mL for mss_cnc
    inputs.GA_units = 1 # glacier algae
    inputs.SA_units = 1 # snow algae
    
    # determine C_factor (can be None or a number)
    # this is the concentrating factor that accounts for
    # resolution difference in field samples and model layers
    inputs.Cfactor_GA = 25
    inputs.Cfactor_SA = 25
    
    # Set names of files containing the optical properties of these LAPs:
    inputs.FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc' # uncoated black carbon (Bohren and Huffman, 1983)
    inputs.FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc' # coated black carbon (Bohren and Huffman, 1983)
    inputs.FILE_brwnC1 = 'brC_Kirch_BCsd.nc' # uncoated brown carbon (Kirchstetter et al. (2004).)
    inputs.FILE_brwnC2 = 'brC_Kirch_BCsd_slfcot.nc' # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.FILE_dust1  = 'dust_balkanski_central_size1.nc' # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    inputs.FILE_dust2  = 'dust_balkanski_central_size2.nc' # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    inputs.FILE_dust3  = 'dust_balkanski_central_size3.nc' # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    inputs.FILE_dust4  = 'dust_balkanski_central_size4.nc' # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    inputs.FILE_dust5 = 'dust_balkanski_central_size5.nc' # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    inputs.FILE_ash1  = 'volc_ash_eyja_central_size1.nc' # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    inputs.FILE_ash2 = 'volc_ash_eyja_central_size2.nc' # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    inputs.FILE_ash3 = 'volc_ash_eyja_central_size3.nc' # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    inputs.FILE_ash4 = 'volc_ash_eyja_central_size4.nc' # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    inputs.FILE_ash5 = 'volc_ash_eyja_central_size5.nc' # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    inputs.FILE_ash_st_helens = 'volc_ash_mtsthelens_20081011.nc' # ashes from Mount Saint Helen's
    inputs.FILE_Skiles_dust1 = 'dust_skiles_size1.nc' # Colorado dust size 1 (Skiles et al 2017)
    inputs.FILE_Skiles_dust2 = 'dust_skiles_size2.nc' # Colorado dust size 2 (Skiles et al 2017)
    inputs.FILE_Skiles_dust3 = 'dust_skiles_size3.nc' # Colorado dust size 3 (Skiles et al 2017)
    inputs.FILE_Skiles_dust4 = 'dust_skiles_size4.nc' # Colorado dust size 4 (Skiles et al 2017)
    inputs.FILE_Skiles_dust5 = 'dust_skiles_size5.nc' # Colorado dust size 5 (Skiles et al 2017)
    inputs.FILE_GreenlandCentral1 = 'dust_greenland_central_size1.nc' # Greenland Central dust size 1 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral2 = 'dust_greenland_central_size2.nc' # Greenland Central dust size 2 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral3 = 'dust_greenland_central_size3.nc' # Greenland Central dust size 3 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral4 = 'dust_greenland_central_size4.nc' # Greenland Central dust size 4 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral5  = 'dust_greenland_central_size5.nc' # Greenland Central dust size 5 (Polashenski et al 2015)
    inputs.FILE_Cook_Greenland_dust_L = 'dust_greenland_Cook_LOW_20190911.nc' # GRIS dust (Cook et al. 2019 "LOW")
    inputs.FILE_Cook_Greenland_dust_C = 'dust_greenland_Cook_CENTRAL_20190911.nc' # GRIS dust 1 (Cook et al. 2019 "mean")
    inputs.FILE_Cook_Greenland_dust_H = 'dust_greenland_Cook_HIGH_20190911.nc' # GRIS dust 1 (Cook et al. 2019 "HIGH")
    inputs.FILE_snw_alg  = 'SA_Halbach2021_r12.4.nc' # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
    inputs.FILE_glacier_algae = 'GA_Chevrollier2021_r5_L17.5.nc' # glacier algae in cells/ml or ppb depending on GA_units (Cook et al. 2020)
    
    
    # Indicate mass mixing ratios scenarios for each impurity (units: ng(species)/g(ice), or ppb)
    # glacier algae in cells/mL if GA_units == 1, ppb if GA_units == 0.
    # The script will loop over the different mixing scenarios
    
    
    inputs.mss_cnc_soot1 = [0]*len(inputs.dz)    
    inputs.mss_cnc_soot2 = [0]*len(inputs.dz)    
    inputs.mss_cnc_brwnC1 = [0]*len(inputs.dz)   
    inputs.mss_cnc_brwnC2 = [0]*len(inputs.dz)   
    inputs.mss_cnc_dust1 = [0]*len(inputs.dz)    
    inputs.mss_cnc_dust2 = [0]*len(inputs.dz)    
    inputs.mss_cnc_dust3 = [0]*len(inputs.dz)    
    inputs.mss_cnc_dust4 = [0]*len(inputs.dz)    
    inputs.mss_cnc_dust5 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash1 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash2 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash3 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash4 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash5 = [0]*len(inputs.dz)    
    inputs.mss_cnc_ash_st_helens = [0]*len(inputs.dz)   
    inputs.mss_cnc_Skiles_dust1 = [0]*len(inputs.dz)   
    inputs.mss_cnc_Skiles_dust2 = [0]*len(inputs.dz)    
    inputs.mss_cnc_Skiles_dust3 = [0]*len(inputs.dz)    
    inputs.mss_cnc_Skiles_dust4 = [0]*len(inputs.dz)  
    inputs.mss_cnc_Skiles_dust5 = [0]*len(inputs.dz)  
    inputs.mss_cnc_GreenlandCentral1 = [0]*len(inputs.dz) 
    inputs.mss_cnc_GreenlandCentral2 = [0]*len(inputs.dz) 
    inputs.mss_cnc_GreenlandCentral3 = [0]*len(inputs.dz) 
    inputs.mss_cnc_GreenlandCentral4 = [0]*len(inputs.dz) 
    inputs.mss_cnc_GreenlandCentral5 = [0]*len(inputs.dz) 
    inputs.mss_cnc_Cook_Greenland_dust_L = [0]*len(inputs.dz) 
    inputs.mss_cnc_Cook_Greenland_dust_C = [0]*len(inputs.dz) 
    inputs.mss_cnc_Cook_Greenland_dust_H = [0]*len(inputs.dz) 
    inputs.mss_cnc_snw_alg = [0,0]    
    inputs.mss_cnc_glacier_algae = [input_concentration,0]   
    
    if write_config_to_textfile:
        omitted_fields = ['tau', 'g', 'SSA', 'mu_not', 'nbr_wvl', 'wvl', 'Fs', 'Fd', 'L_snw', 'flx_slr']
    
        with open('./model_config.txt', 'w') as f:
            for i in inputs._fields:
                if i not in omitted_fields:
                    f.write(i)
                    f.write(':  ')
                    f.write(str(getattr(inputs,str(i))))
                    f.write('\n')
    
    ##########################################################################
    ################## CALL FUNCTIONS AND PLOT OUTPUTS #######################
    ##########################################################################

    #########################################
    ## IF NO INPUT ERRORS --> FUNCTION CALLS
    #########################################
    
    
    outputs = snicar_feeder(inputs)
    albedo = outputs.albedo 

    
    return albedo  

def gen_LUT(densities, grain_sizes, dz, layer_types, sza, path):
    LUT=pd.DataFrame()
    for layer_type in layer_types:
        for density in densities:
            for radius in grain_sizes:
                for depth in dz:
                    for angle in sza:
                        albedo = call_SNICAR(depth, 0, density, layer_type, radius, angle) # albedo between 0.2 and 5µm with 10nm step
                        LUT[f'density={density} kg m3, depth={depth[1]} m, bbl_size={radius} um, sza={angle}, layer_type={layer_type}']=albedo
    LUT.to_pickle(path)


input_concentration = 0 #26x45x24x3 
densities= [300, 320, 340]#, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580,600,
            #620, 640, 660, 680, 700, 720, 740, 760, 780, 800] #density of first layer
grain_sizes = [1000,1500]# ,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,
               #8500, 9000,9500, 10000, 10500, 11000, 11500, 12000, 12500,13000, 13500,14000,
               #14500, 15000,15500,16000, 16500, 17000,17500, 18000, 18500,19000,19500, 20000]
dz =[[0.001, 0.01],[0.001, 0.015]]#, [0.001,0.02],[0.001,0.025],[0.001,0.03],[0.001,0.035], 
     #[0.001,0.04],[0.001,0.045], [0.001,0.05],[0.001,0.055],[0.001,0.06],[0.001,0.065], 
     #[0.001,0.07],[0.001,0.075],[0.001,0.08],[0.001,0.085], [0.001,0.09],[0.001,0.095], 
     #[0.001,0.1], [0.001,0.105],[0.001,0.11],[0.001,0.115],[0.001,0.120],[0.001,0.125]] 
layer_types = [[0,0]] #has to be the same length than the depths in dz
sza=[35,40,45]
gen_LUT(densities, grain_sizes, dz, layer_types, sza, '/Users/au660413/Desktop/LUT_inversemodelling_051121_00.pkl')
