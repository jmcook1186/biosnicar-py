#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joe, lou

Driver for the bio-optical model developed by Cook et al. 2017,2020 to
calculate algae ACS (absorption cross section). The model workflow consists
in three functions detailed below in case the user needs to use the functions
separately. If the user wants to run the full model, the parameters to document
are summarized at the beginning of the code.

######################################################################################
Model workflow
######################################################################################

    1. bioptical_calculations() : Calculation of absorption cross section (ACS)
        in m2/cell or m2/um3 and refractive index (n,k; unitless) in the
        spectral range of interest, both at 1nm and 10nm resolution
        (used in BioSNICAR)
    2. ssp_calculations() : Calculation of single scattering properties
        (g, ssa; unitless) using Mie or geometric optics theory in the
        spectral range of interest at the resolution of BioSNICAR (10nm)
    3. save_netcdf() : Storage of the data and associated metadata in a netcdf
        file directly usable in BioSNICAR

######################################################################################
Inputs of the different functions
######################################################################################

    1. wvl = (numpy array) spectral range of interest (in µm, 1nm step)
       ACS_calculated = (boolan) True if the ACS is calculated from pigment
                        data below, False if loaded
       ACS_file = (string - if ACS_calculated false) directory to the ACS
                    file (m2/cell)
       pigment_dir = (string) directory to file containing pigment mass
                    absorption coefficients - must be csv file with size
                    and resolution of wvl, and coefficients in m2/mg
       pigments_data = (string) dictionary with pigment file names and
                        associated concentrations (ng/cell or mg/µm3)
       biovolume = (boolean) True if concentrations of pigments in mg/µm3
                    (ACS generated in m2/um3), False if given in ng/cell
                    (ACS generated in m2/cell)
       density =  (int - used if biovolume = True) density of wet biomass
                   (kg/µm3 - 1160*10**(-18) dflt value from green snw algae)
       xw = (float - not used atm) water fraction of a cell
                    (0.7-0.8 by default)
       cell_volume = (int - used if biovolume = False) volume of an algae cell
                     reference values from Halbach et al. 2021:
                     1500 for GA, 9600 for SA
       n_algae = (numpy array) real part of cellular refractive index
                  in the spectral range of wvl (constant 1.4 by default,
                  Dauchet et al. 2015)
       k_water = (numpy array) imaginary part of the refractive index of water
                  in the spectral range of wvl
       packaging_correction = (boolean - only if ACS_calculated is True)
                            if True, absorption coefficients retrieved
                            from summing pigments are linearly corrected
                            (Williamson et al. 2020)
       smooth = (boolean) if True,  apply optional smoothing filter
       window_size = (int) size of window of smoothing filter (dflt value: 25)
       poly_order = (int) polynomial order of smoothing filter (dflt value: 3)
       smoothStart = (int) start of smoothing filter (default value: 44)
       smoothStop = (int) stop of smoothing filter (default value: 100)
       plot_n_k_ACS_cell = (boolean) if True, plot with n,k and ACS printed
       plot_pigment_ACS = (boolean) if True, plot with pigment ACSs printed
       savefiles_n_k_ACS_cell = (boolean) if True, files with k,n and ACS
                                saved in the directory savepath_n_k_ACS_plots
       saveplots_n_k_ACS = (boolean) if True, plots saved in the directory
                            savepath_n_k_ACS_plots
       savefilename_n_k_ACS_cell = (string) name of the file containing n,
                                k and ACS if savefiles toggled on
       savepath_n_k_ACS_plots = (boolean) directory for saving data if
                                savefiles or saveplots toggled on

     2. GO = (boolean) if True, uses geometric optics to calculate single
            scattering OPs assuming cell shape = cylinder
        Mie = (boolean) if True, uses Mie theory to calculate single
            scattering OPs assuming cell shape = sphere
        r = (int) radius of sphere (Mie)/cynlinder (GO) representing cell (µm)
        L = (int) depth of the cylinder representing the cell (GO option, µm)
        wvl = (numpy array) same as in 1.
        n_algae = (numpy array) with real part of the RI of a single cell
            in the spectral range of wvl
        k_algae = (numpy array) with imaginary part of RI of a single cell
             in the spectral range of wvl (output of bioptical_calculations())
        plot_OPs = (boolean) if True, print plots with OPs
        savefigs_OPs = (boolean) if True, save plots with OPs
        figname_OPs = (string) figure name
        report_dims = (boolean) if True, cell dimensions printed to console

     3. netcdf_save = (boolean) if True, saves data in a netcdf file
        GO = (boolean) if True, uses geometric optics to calculate single
            scattering OPs assuming cell shape = cylinder
        Mie = (boolean) if True, uses Mie theory to calculate single
            scattering OPs assuming cell shape = sphere
        g = (numpy array) assym parameter retrieved by ssp_calculations()
        ssa = (numpy array) assym parameter retrieved by ssp_calculations()
        ACS = (numpy array) mass absorption coefficient retrieved by
            bioptical_calculations()
        L = (numpy array) used only if GO set to True, depth of the cylinder
            representing the cell (µm)
        r = (numpy array) radius of the sphere/cynlinder representing the
            cell (µm)
        density = (int) density of dry or wet biomass (kg/µm3 - 1400*10**(-18)
                for dry biomass, Dauchet et al. 2015)
        savepath_netcdf = (string) save path directory
        filename_netcdf = (string) name of the file containing the
                        optical properties
        information = (string) containing any additional info for
                    metadata in netcdf (e.g. 'Glacier algae OPs derived
                    from GO calculations with empirical ACS')

######################################################################################
Outputs of the different functions
######################################################################################
    1. numpy arrays for wvl, n, k and CAC in the spectral range 200-5000µm w/
        both 1nm resolution and 10nm (=BioSNICAR) resolution (xxx_rescaled)
        for a given cell;
        dataframe combining variables at BioSNICAR (10nm) resolution;
        dataframe with individual pigment absorption coefficients
    2. numpy arrays with ssa and g for a single cell in the spectral
        range 200-5000µm at BioSNICAR (10nm) resolution
    3. netcdf file with CAC, ssa and g that can be used directly in BioSNICAR
    
    
Note:
In GO mode, the main calculations are based upon the equations of
Diedenhoven et al (2014), who provided a python script as supplementary
material for their paper, available at:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

In Mie mode, the optical properties are calculated using Mie scattering
using Scott Prahl's miepython package https://github.com/scottprahl/miepython.
"""
#%%
from BioOptical_Model.biooptical_Funcs import bioptical_calculations, ssp_calculations, net_cdf_updater
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import xarray as xr
from netCDF4 import Dataset
import collections as clt
from scipy.interpolate import interp1d
import seaborn as sns

#%%

######################################################################################
## INPUTS TO FILL
######################################################################################

######## Set spectral range and base directory pointing to package
dir_base = '/home/joe/Code/BioSNICAR_GO_PY/'
wvl = np.arange(0.2, 4.999, 0.001) #spectral range of interest in µm

######## if ACS is calculated:
ACS_calculated = False
biovolume = False
k_water = np.loadtxt(dir_base + 'Data/pigments/k_ice_480.csv')
packaging_correction = False
pigment_dir =  dir_base + 'Data/pigments/'
pigments_data = {str(pigment_dir + 'alloxanthin.csv'): 0.0,
                 str(pigment_dir + 'antheraxanthin.csv'): 0,
                 str(pigment_dir + 'chl-a.csv'): 3.96e-3,
                 str(pigment_dir + 'chl-b.csv'): 7e-4,
                 str(pigment_dir + 'lutein.csv'): 0,
                 str(pigment_dir + 'neoxanthin.csv'): 0,
                 str(pigment_dir + 'pheophytin.csv'): 0.0,
                 str(pigment_dir + 'Photop_carotenoids.csv'): 0.0,
                 str(pigment_dir + 'Photos_carotenoids.csv'): 6e-3,
                 str(pigment_dir + 'ppg_shifted.csv'): 4.3e-2,
                 str(pigment_dir + 'trans_astaxanthin_ester.csv'): 0.0,
                 str(pigment_dir + 'trans_astaxanthin.csv'): 0,
                 str(pigment_dir + 'violaxanthin.csv'): 0,
                 str(pigment_dir + 'zeaxanthin.csv'): 0,
                }

density = 1160*10**(-18)
xw = 0.8

# Optional smoothing filter for calculated k, ACS
smooth = True
window_size = 25
poly_order = 3
smoothStart = 44
smoothStop = 100

######## if ACS is loaded from a file:
ACS_file = '/home/joe/Code/BioSNICAR_GO_PY/Data/ACS_GA_Halbach2021.csv'
#FIELDWORK1/SA_ACS_Chevrollier_2021.csv'


######## Chose method for calculation of optical properties
GO = False
Mie = True

######## Algae properties for calculations of optical properties
n_algae = 1.4 * np.ones(np.size(wvl)) 
r = 5
L = 15 #5 and 15 Halbach 2021.
cell_volume = 1178# 7900 # Chevrollier et al. 2021: SA 2572, GA 1374; Halbach et al. 2021: 1178 GA, 7900 SA


####### Directories and printing/saving options for calculated k, ACS
plot_n_k_ACS_cell = True
plot_pigment_ACSs = False
savefiles_n_k_ACS_cell = False
saveplots_n_k_ACS = False
savepath_n_k_ACS_plots = dir_base
savefilename_n_k_ACS_cell = ''

######## Directories and printing/saving options for optical properties
plots_OPs = True
savefigs_OPs = False
report_dims = False
savepath_OPs = dir_base
figname_OPs = 'figname'

######## Saving OPs in netcdf
netcdf_save = True
savepath_netcdf = dir_base + 'Data/Mie_files/480band/lap/'
filename_netcdf = 'GA_Halbach2021_oct21'
information = ''

#%%
##############################################################################
## CALCULATIONS OF ABSORPTION PROPERTIES
##############################################################################

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR,\
ACS, ACS_rescaled_BioSNICAR, n_k_ACS_rescaled_BioSNICAR,\
abs_coeff_pigm_DataFrame = bioptical_calculations(ACS_calculated, ACS_file,
                                                  biovolume, density, xw,
                                                  cell_volume, wvl,
                                                  packaging_correction,
                                                  pigment_dir, pigments_data,
                                                  n_algae, k_water, smooth, 
                                                  window_size, poly_order,
                                                  smoothStart, smoothStop,
                                                  plot_n_k_ACS_cell,
                                                  plot_pigment_ACSs,
                                                  savefiles_n_k_ACS_cell,
                                                  savefilename_n_k_ACS_cell,
                                                  saveplots_n_k_ACS,
                                                  savepath_n_k_ACS_plots)

    
#%%
##############################################################################
## CALCULATIONS OF SCATTERING PROPERTIES
##############################################################################

assym, ss_alb = ssp_calculations(GO, Mie, savepath_OPs, r, L, 
                                 wvl_rescaled_BioSNICAR,
                                 n_rescaled_BioSNICAR,
                                 k_rescaled_BioSNICAR,
                                 plots_OPs, savefigs_OPs,
                                 figname_OPs, report_dims)

#%% 
##############################################################################
## SAVING DATA IN NETCDF
##############################################################################

if netcdf_save:
    net_cdf_updater(GO, Mie, savepath_netcdf, filename_netcdf,
                    wvl_rescaled_BioSNICAR, assym, ss_alb,
                    ACS_rescaled_BioSNICAR, L, r, density, information)
