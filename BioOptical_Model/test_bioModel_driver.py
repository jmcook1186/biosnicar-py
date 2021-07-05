#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TO DO : reduce to column 79
"""
Created on Sat Apr 17 15:45:55 2021

@author: joe, lou

Driver for the bio-optical model developed by Cook et al. 2020 adapted to calculate CAC (cell absorption coefficient) as in Halbach et al. 2021.
The model workflow consists in three functions detailed below in case the user needs to use the functions separately. 
If the user wants to run the full model, the parameters to document are summarized at the beginning of the code.

######################################################################################
Model workflow
######################################################################################

    1. bioptical_calculations() : Calculation of cell absorption coefficient (CAC) and refractive index (n,k) in the spectral range of interest
    2. ssp_calculations() : Calculation of single scattering properties (g, ssa) using Mie or geometric optics theory in the spectral range of interest
    3. save_netcdf() : Storage of the data and associated metadata in a netcdf file directly usable in BioSNICAR

######################################################################################
Inputs of the different functions
######################################################################################

    1. wvl = (numpy array) spectral range of interest (in µm)
       MAC_calculated = (boolan) True if the MAC is calculated from pigment data below, False if loaded
       MAC_path = (string) directory to the MAC file if loaded, units in m2/cell 
       pigment_dir = (string) directory to pigment mass absorption coefficients (must be csv files with wvl size and resolution, in m2/mg)
       pigments_data = (string) dictionary with pigment file names and associated concentrations (ng/cell or mg/µm3)
       biovolume = (boolean) True if concentrations of pigments given in mg/µm3, False if given in ng/cell
       density =  (int) density of wet biomass (kg/µm3 - 1160*10**(-18) from lab trials on snow algae), only for biovolume=True (!)
       xw = (float) water fraction of a cell (0.7-0.8 by defaultp)
       cell_volume = (int) volume of an algae cell, only for biovolume=False (!)
       n_algae = (numpy array) real part of the refractive index of a single cell in the spectral range of wvl (constant 1.4 by default, Dauchet et al. 2015)
       k_water = (numpy array) imaginary part of the refractive index of water in the spectral range of wvl
       packaging_correction = (boolean) if True, absorption coefficients are linearly corrected (Williamson et al. 2020?)
       smooth = (boolean) if True,  apply optional smoothing filter 
       window_size = (int) size of window of smoothing filter
       poly_order = (int) order of polynomial of smoothing filter
       smoothStart = (int) start of smoothing filter
       smoothStop = (int) stop of smoothing filter
       plot_n_k_MAC_cell = (boolean) if True, returns plot with n,k and MAC
       plot_pigment_MACs = (boolean) if True, returns plot with MAC of different pigments
       savefiles_n_k_MAC_cell = (boolean) if True, save n, k and MAC in the pigment directory as csv files
       savefilename_n_k_MAC_cell = (string) name of the file containing n, k and MAC if savefiles toggled on
       saveplots_n_k_MAC = (boolean) if True, save plots in the pigment directory
       savepath_n_k_MAC = (boolean) directory for saving data if savefiles or saveplots toggled on

       
     2. GO = (boolean) if True, uses geometric optics to calculate single scattering OPs assuming cell shape = cylinder
        Mie = (boolean) if True, uses Mie theory to calculate single scattering OPs assuming cell shape = sphere
        r = (int) radius of the sphere/cynlinder representing the cell (µm)
        L = (int) depth of the cylinder representing the cell (GO option, µm)
        wvl = (numpy array) wavelengths in spectral range used in bioptical calculations (µm)
        n_algae = (numpy array) with real part of the RI of a single cell in the spectral range of wvl (constant 1.4 by default, Dauchet et al. 2015)
        k_algae = (numpy array) with imaginary part of RI of a single cel in the spectral range of wvl (output of bioptical_calculations())
        plot_OPs = (boolean) if True, returns plots with OPs
        savefigs_OPs = (boolean) if True, save plots with OPs
        figname_OPs = (string) figure name
        report_dims = (boolean) if True, cell dimensions printed to console
        
     3. netcdf_save = (boolean) if True, saves data in a netcdf file  
        g = (numpy array) assym parameter retrieved by ssp_calculations()
        ssa = (numpy array) assym parameter retrieved by ssp_calculations()
        MAC = (numpy array) mass absorption coefficient retrieved by bioptical_calculations()
        L = (numpy array) used only if GO set to True, depth of the cylinder representing the cell (µm)
        r = (numpy array) radius of the sphere/cynlinder representing the cell (µm)
        density = (int) density of dry or wet biomass (kg/µm3 - 1400*10**(-18) for dry biomass, Dauchet et al. 2015)
        savepath_netcdf = (string) save path directory 
        filename_netcdf = (string) name of the file containing the optical properties
        information = (string) containing any additional info for metadata in netcdf (e.g. 'Optical properties derived from geometrical optics calculations with empirically derived MAC')

######################################################################################
Outputs of the different functions
######################################################################################
    1. numpy arrays for wvl, n, k and CAC in the spectral range 200-5000µm at 1 and 10nm resolution for a given cell; 
        dataframe combining variables at BioSNICAR (10nm) resolution; 
        dataframe with individual pigment absorption coefficients
    2. numpy arrays with ssa and g for a single cell in the spectral range 200-5000µm at BioSNICAR (10nm) resolution
    3. netcdf file with CAC, ssa and g that can be used directly in BioSNICAR
    
    
Note: 
In GO mode, the main calculations are based upon the equations of Diedenhoven et al (2014),
who provided a python script as supplementary material for their paper, available at:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

In Mie mode, the optical properties are calculated using Mie scattering using Scott Prahl's miepython
package https://github.com/scottprahl/miepython. 
"""
#%%
from test_bioModel import bioptical_calculations, ssp_calculations, net_cdf_updater
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

######## Set spectral range
wvl = np.arange(0.2,4.999,0.001) #spectral range of interest in µm

######## Enter algae properties 
MAC_calculated = True
MAC_path = '/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/MAC_HE2C.csv'
biovolume = False #False if MAC in m2/cell, True if MAC in m2/µm3 (depends on the pigment concentrations)
density = 1160*10**(-18) #density of wet biomass determined on snow algae (kg/µm3 - 1400*10**(-18) for dry biomass, Dauchet et al. 2015)
pigment_dir = '/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/' #directory with pigment mass absorption coefficients in the range of wvl, in m2/kg

#### pigment data: name of the file in the above directory + associated concentration in ng/cell or mg/µm3
pigments_data = ''

xw = 0.8 #water fraction of a cell (0.7-0.8 by default)
r = 9 #radius of the cell in µm0
L = 18 #length of the cell in µm
cell_volume = 1900 # volume of single cell in µm3 used for biovolume=True, ex 1500 for MB/AN, 9600 for CN
n_algae = 1.4*np.ones(np.size(wvl)) #np array of real part of the refractive index of a single cell in the spectral range of interest (constant 1.4 by default, Dauchet et al. 2015)
k_water = np.loadtxt('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/k_ice_480.csv')  # np array with imaginary part of the refractive index of water in the spectral range of interest
packaging_correction = False #if true, absorption coefficients are linearly corrected

####### Optional smoothing filter for calculated k,MAC
smooth = True #if true,  apply optional smoothing filter with
window_size = 25 #size of window of smoothing filter ()
poly_order = 3 #order of polynomial of smoothing filter
smoothStart = 44 #start of smoothing filter (44 default)
smoothStop = 100 #stop of smoothing filter (100 default)

####### Directories and file/plot saving for calculated k, MAC
savepath_n_k_MAC = '/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/' #directory for saving ,k and MAC if savefiles or saveplots toggled on
plot_n_k_MAC_cell = True #if true, returns plot with n,k and MAC of single cell
plot_pigment_MACs = False #if true, returns plot with MAC of different pigments
savefiles_n_k_MAC_cell = False #if true, save n, k and MAC of single cell in the savepath as csv files
savefilename_n_k_MAC_cell  = 'avg_snow_algae_wout_esters' #name of the file containing n,k and MAC if savefiles toggled on
saveplots_n_k_MAC = False #if true, save plots in the savepath

######## Chose method for calculation of optical properties
GO = True #if True, uses geometric optics to calculate single scattering OPs assuming cell shape = cylinder
Mie = False #if True, uses Mie theory to calculate single scattering OPs assuming cell shape = sphere

######## Plotting/printing options for optical properties
savepath_OPs=pigment_dir
plots_OPs = False #True to plot optical properties
savefigs_OPs = False #True to save the figures in the savepath
figname_OPs = '' #string with figure name
report_dims = False #Boolean toggling cell dimensions printing to console

######## Saving OPs in netcdf
netcdf_save = False
savepath_netcdf = '/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/Mie_files/480band/lap/' # directory where the netcdf file should be saved
filename_netcdf = 'community_BR' # name of the netcdf file
information = 'from Halbach et al. 2021' #string containing information as metadata in the netcdf file

#%%
######################################################################################
## CALCULATIONS OF ABSORPTION PROPERTIES
######################################################################################
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
    
#%%
######################################################################################
## CALCULATIONS OF SCATTERING PROPERTIES
######################################################################################

assym, ss_alb = ssp_calculations(GO, Mie, savepath_OPs, r, L, wvl_rescaled_BioSNICAR, n_rescaled_BioSNICAR, k_rescaled_BioSNICAR, plots_OPs, savefigs_OPs, figname_OPs, report_dims)

#%% 
######################################################################################
## SAVING DATA IN NETCDF
######################################################################################

if netcdf_save:
    net_cdf_updater(savepath_netcdf, filename_netcdf, wvl_rescaled_BioSNICAR, assym, ss_alb, MAC_rescaled_BioSNICAR, L, r, density, information)



#%%
################################################################################################################################################
################################################################################################################################################
### TESTINGS
################################################################################################################################################
################################################################################################################################################



#%% RAW EXTRACT ABS COEFF CALCULATION
#### Generation of ppg abs coeff from Laura's data with raw extracts 

df_as=pd.read_csv('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/ppg_raw_extract_all_stations.csv',index_col='wvl')
df_as=df_as.loc[:,~df_as.columns.str.contains('MIT1')] #remove snow sites
df_as=df_as.loc[:,~df_as.columns.str.contains('MIT2')] #remove snow sites
wvl_interp = np.arange(200,5000,1) #spectral range of interest 
df_interpolated=pd.DataFrame(index=wvl_interp)
for col in df_as.columns:
    df_as[col]=2.303*df_as[col]
    f_interp = interp1d(df_as.index,df_as[col], kind='cubic',fill_value='extrapolate')
    df_interpolated[col]=f_interp(wvl_interp)
    df_interpolated[col]=df_interpolated[col]-df_interpolated[col][750] #baselined at 750nm
df_interpolated.loc[750:, :]=0
df_interpolated.loc[:280, :]=0
coeff_average=df_interpolated.mean(axis=1)
coeff_average=coeff_average-coeff_average[800]
test=df_interpolated.plot(xlim=(200,800), ylim=(0, 0.05), title='ppg abs coeff derived from raw extracts')
test.set_ylabel("abs (m$^2$ mg$^{-1}$)")
#test.figure.savefig('/Users/au660413/Desktop/ppg_abs.png', format='png')
#coeff_average.to_csv('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/ppg_raw_extract_trial.csv',header=None,index=False)

#%% Calculation MAC Laura data for all stations 

MAC_df=pd.DataFrame(index=wvl*1000)
names=['BR', 'HE1','HE2','MIT6','MIT1','MIT2','MIT3', 'MIT4', 'MIT5']
pigment_data_list=[pigments_data_BR,pigments_data_HE1, pigments_data_HE2,
                   pigments_data_glacier_algae_MIT6_community,
                   pigments_data_MIT1,pigments_data_MIT2,pigments_data_MIT3,
                   pigments_data_MIT4,pigments_data_MIT5]


r_list=9*np.ones(len(pigment_data_list))
L_list=18*np.ones(len(pigment_data_list))
#savepath_netcdf=/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/

for pigments_data,r,L,name in zip(pigment_data_list, r_list,L_list, names):
    wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR, data, abs_coeff \
        = bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
    MAC_df[name]=MAC #m2/cell
    assym, ss_alb = ssp_calculations(GO, Mie, savepath_OPs, r, L, wvl_rescaled_BioSNICAR, n_rescaled_BioSNICAR, k_rescaled_BioSNICAR, plots_OPs, savefigs_OPs, figname_OPs, report_dims)
    net_cdf_updater(savepath_netcdf, name, wvl_rescaled_BioSNICAR, assym, ss_alb, MAC_rescaled_BioSNICAR, L, r, density, information)
        
#test=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/BioOptical_Model/algae_geom_4_40.nc')
#MAAC_plot=MAC_df.loc[:,~MAC_df.columns.str.contains('shifted')] #select only non shifted
#MAC_df['Cook 2020']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(1*10**12))])
#ax = MAC_df.plot(use_index=True, xlim=(300,800), title='Glacier algae cellular MAC - Sermilik 2019')
#plt.legend(fontsize='large', ncol=2,handleheight=2.4, labelspacing=0.05)
#ax.set_xlabel("wvl (nm)")
#ax.set_ylabel("absorption coefficient (m$^2$ cell$^{-1}$)")
#ax.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/community_reconstruction.pdf')
MAC_df.to_csv('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/CAC_data.csv')

#%% Calculation MAC Laura data min mean max 

MAC_df_standard=pd.DataFrame(index=wvl[::10]*1000)
MAC_df_extract=pd.DataFrame(index=wvl[::10]*1000)
#r_list=[9.1,9.4,10.1]
#L_list=[18.4,17,23]

names=['min_MIT5a','average','max_MIT3b']
pigment_data_list=[pigments_data_glacier_algae_BR1_shifted, pigments_data_glacier_algae_BR2_shifted, pigments_data_glacier_algae_BR3_shifted,\
                   pigments_data_glacier_algae_HE2a_shifted,pigments_data_glacier_algae_HE2b_shifted,pigments_data_glacier_algae_HE2c_shifted,\
                   pigments_data_glacier_algae_MIT6_shifted,\
                   pigments_data_glacier_algae_MIT3a_shifted,pigments_data_glacier_algae_MIT3b_shifted,pigments_data_glacier_algae_MIT3c_shifted,\
                   pigments_data_glacier_algae_MIT4a_shifted,pigments_data_glacier_algae_MIT4b_shifted,pigments_data_glacier_algae_MIT4c_shifted,\
                   pigments_data_glacier_algae_MIT5a_shifted, pigments_data_glacier_algae_MIT5b_shifted, pigments_data_glacier_algae_MIT5c_shifted]

pigment_data_min_mean_max=[pigments_data_glacier_algae_MIT5a_shifted,average,pigments_data_glacier_algae_MIT3b_shifted]

for pigments_data,r,L in zip(pigment_data_min_mean_max_standard, r_list,L_list):
    wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR, data, abs_coeff \
        = bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
    MAC_df_standard['']=MAC_rescaled_BioSNICAR #m2/cell

pigment_data_list=[pigments_data_glacier_algae_BR1_extract, pigments_data_glacier_algae_BR2_extract, pigments_data_glacier_algae_BR3_extract,\
                   pigments_data_glacier_algae_HE2a_extract,pigments_data_glacier_algae_HE2b_extract,pigments_data_glacier_algae_HE2c_extract,\
                   pigments_data_glacier_algae_MIT6_extract,\
                   pigments_data_glacier_algae_MIT3a_extract,pigments_data_glacier_algae_MIT3b_extract,pigments_data_glacier_algae_MIT3c_extract,\
                   pigments_data_glacier_algae_MIT4a_extract,pigments_data_glacier_algae_MIT4b_extract,pigments_data_glacier_algae_MIT4c_extract,\
                   pigments_data_glacier_algae_MIT5a_extract, pigments_data_glacier_algae_MIT5b_extract, pigments_data_glacier_algae_MIT5c_extract]

#Calculate average accross stations
total=clt.Counter({})
average=clt.Counter({})
nb_stations=len(pigment_data_list)
for dictionary in pigment_data_list:
    total=total+clt.Counter(dictionary)
for key in total:
    average[key]=total[key]/nb_stations
    
pigment_data_min_mean_max=[pigments_data_glacier_algae_MIT5a_shifted,average,pigments_data_glacier_algae_MIT3b_shifted]
 
for pigments_data,r,L in zip(pigment_data_min_mean_max_extract, r_list,L_list):
    wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR, data, abs_coeff \
        = bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
    MAC_df_extract['']=MAC_rescaled_BioSNICAR #m2/cell
    
MAC_df_standard.columns=names
MAC_df_extract.columns=names
#bigdata=pd.concat([MAC_df_standard, MAC_df_extract], axis=1)
#bigdata["wvl"]=MAC_df_extract.index
#bigdata.columns=['extract', 'extract', 'extract', 'reconstruct', 'reconstruct','reconstruct', 'wvl']
#tryy=pd.melt(MAC_df, id_vars="wvl")
g=sns.lineplot(data=tryy, x="wvl", y="value", hue="variable")
g.set(xlim=(300, 800), ylabel="MAC (m$^2$ cell$^{-1}$)", xlabel="wavelength (nm)", title='Cellular MAC - Sermilik 2019')
#g.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/extract_vs_reconstruction.eps', format='eps')

#%%
test=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/BioOptical_Model/algae_geom_4_40.nc')
MAC_plot=MAC_df.loc[:,~MAC_df.columns.str.contains('shifted')] #select only non shifted
MAC_plot['Cook 2020']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(1*10**12))])
ax = MAC_plot.plot(use_index=True, xlim=(300,800), title='Glacier algae cellular MAC - Sermilik 2019')
ax.set_xlabel("wvl (nm)")
ax.set_ylabel("MAC (m$^2$ cell$^{-1}$)")
ax.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/min_mean_max_MAC_with_extract.eps', format='eps')

#%% Comparison IS data with Cook 2020

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_shift, data, abs_coeff \
        = bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_HE2c_shifted,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_extract, data, abs_coeff \
        = bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_HE2c_extract,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)

MAC_plot=pd.DataFrame(index=wvl[::10]*1000)
IS_data=pd.read_csv('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/IS_data.csv')
#HE2C = np.array(IS_data['HE2C']).flatten()
#MIT3A = np.array(IS_data['MIT3A']).flatten()
#MIT3C = np.array(IS_data['MIT3C']).flatten()
#MIT1A = np.array(IS_data['MIT1A']).flatten()
#MIT1B = np.array(IS_data['MIT1B']).flatten()
#HE2A = np.array(IS_data['HE2A']).flatten()
#MIT5C = np.array(IS_data['MIT5C']).flatten()
#MIT5A = np.array(IS_data['MIT5A']).flatten()
#MIT5B = np.array(IS_data['MIT5B']).flatten()
#test=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/BioOptical_Model/algae_geom_4_40.nc')
#MAC_plot['Cook 2020']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(1*10**12))])
#MAC_plot['IS_HE2C']=HE2C[::10]
#MAC_plot['IS_MIT3A']=MIT3A[::10]
#MAC_plot['IS_MIT3C']=MIT3C[::10]
#MAC_plot['IS_MIT1A']=MIT1A[::10]
#MAC_plot['IS_MIT1B']=MIT1B[::10]
#MAC_plot['IS_HE2A']=HE2A[::10]
#AC_plot['IS_MIT5C']=MIT5C[::10]
#MAC_plot['IS_MIT5A']=MIT5A[::10]
#MAC_plot['IS_MIT5B']=MIT5B[::10]
IS_average=np.array(IS_data.mean(axis=1)).flatten()
MAC_plot['IS_average']=IS_average[::10]
MAC_plot['Extract_mean']=np.array(MAC_df.mean(axis=1)).flatten()
#MAC_plot['Standard_mean']=MAC_rescaled_BioSNICAR_shift
ax = MAC_plot.plot(use_index=True, xlim=(300,800), title='Glacier algae cellular MAC - Sermilik 2019')
ax.set_xlabel("wvl (nm)")
ax.set_ylabel("MAC (m$^2$ cell$^{-1}$)")
#ax.figure.savefig('/Users/au660413/Documents/PhD/Data/Laura_pigments/IS_pigment_comparison_MIT3c.eps', format='eps')


#%%
# Comparison MAC using Chris's data

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_all_pigs, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_allpig,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_phenols_original, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_phenols_original,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
        
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_phenols_corrected, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                        packaging_correction,pigment_dir, pigments_data_phenols_corrected,\
                        n_algae, k_water, smooth, window_size,\
                        poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                        plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                        savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_phenols_invivocorrected, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_phenols_invivocorrected,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_ppg, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_ppg,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                       savepath_n_k_MAC)

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_ppg_shifted, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_ppg_shifted,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                       savepath_n_k_MAC)

wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_ppg_extract, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, average,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                       savepath_n_k_MAC)
    

MAC_df=pd.DataFrame(index=wvl_rescaled_BioSNICAR*1000)
#MAC_df['Chris_MAC_calculated_allpig_phenols']=MAC_rescaled_BioSNICAR_all_pigs #from m2/cell to m2/kg dry weight with cell=1.1ng
phenol_correction=(np.array(pd.read_csv('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/pigments_trials/phenol_mac_correction.csv'))).flatten() 
#MAC_df['phenol_invivo_corrected']=MAC_rescaled_BioSNICAR_phenols_invivocorrected
#MAC_df['Chris_MAC_calculated_ppg']=MAC_rescaled_BioSNICAR_ppg #from m2/cell to m2/kg dry weight with cell=1.1ng
#MAC_df['ppg_from_Laura']=MAC_rescaled_BioSNICAR_ppg_shifted #from m2/cell to m2/kg dry weight with cell=1.1ng
test=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/BioOptical_Model/algae_geom_4_40.nc')
#MAC_df['Williamson et al. 2020 - uncorrected']=MAC_rescaled_BioSNICAR_phenols_original
#MAC_df['Williamson et al. 2020 - corrected']=MAC_rescaled_BioSNICAR_phenols_corrected
#MAC_df['Cook 2020 with 1 cell = 1ng']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(1*10**(12)))])
#MAC_df['Cook 2020 with 1 cell = 0.84ng']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(0.84*10**(12)))])
MAC_df['this study']=MAC_rescaled_BioSNICAR_ppg_extract #from m2/cell to m2/kg dry weight with cell=1.1ng
MAC_algae_github=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/Mie_files/480band/lap/Glacier_Algae_480.nc')
MAC_df['C. Williamson unpubl. data with 1 cell = 1ng']=(MAC_algae_github['ext_cff_mss'])/(1*10**(12))
#MAC_df.columns=['Cook 2020 with 1 cell = 1ng', 'Williamson 2020 - corrected', 'this study']

### PLOT comparison.Chris
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(MAC_df.index, MAC_df['Williamson et al. 2020 - uncorrected'], 'g-')
# ax1.plot(MAC_df.index, MAC_df['Williamson et al. 2020 - corrected'], 'g--')
# ax2.plot(MAC_df.index, MAC_df['this study'], 'b-')
# ax2.set_ylim((0,25e-10))
# ax1.set_ylim((0,2.5e-8))
# ax1.set_xlim((280,800))
# ax2.set_xlim((280,800))
# ax1.set_xlabel("wavelength (nm)")
# ax1.set_ylabel("MAC (m$^2$ cell$^{-1}$)", color="green")
# ax2.set_ylabel("MAC (m$^2$ cell$^{-1}$)", color="blue")
# ax1.legend(['Williamson et al. 2020 - uncorrected', 'Williamson 2020 - corrected'])
# ax2.legend(['this study'], loc='upper right', bbox_to_anchor=(0.737, 0.9))
# #ax2.legend(['','','this study'])
# ax2.grid(None)

### PLOT comparison.Cook2020
# ax1=MAC_df.plot(xlim=(300,800))
# ax1.plot(MAC_df.index, MAC_df['Williamson et al. 2020 - corrected'], 'g--')
# ax1.set_xlabel("wavelength (nm)")
# ax1.set_ylabel("MAC (m$^2$ cell$^{-1}$)")

### PLOT comparison.unpubldata
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(MAC_df.index, MAC_df['C. Williamson unpubl. data with 1 cell = 1ng'], 'g-')
ax2.plot(MAC_df.index, MAC_df['this study'], 'b-')
ax2.set_ylim((0,1e-9))
ax1.set_ylim((0,1e-8))
ax1.set_xlim((280,800))
ax2.set_xlim((280,800))
ax1.set_xlabel("wavelength (nm)")
ax1.set_ylabel("MAC (m$^2$ cell$^{-1}$)", color="green")
ax2.set_ylabel("MAC (m$^2$ cell$^{-1}$)", color="blue")
ax1.legend(['C. Williamson unpubl. data with 1 cell = 1ng'])
ax2.legend(['this study'], loc='upper right', bbox_to_anchor=(0.66, 0.95))
#ax2.legend(['','','this study'])
ax2.grid(None)

ax1.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/MAC_comparison_CWunpubl.eps', format='eps')

#%%
abs_cff_pigments=pd.DataFrame(index=wvl*1000)

pigments_data = {str(pigment_dir+'ppg.csv'):1.78e-4,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.97e-2,
                  str(pigment_dir+'phenols_original.csv'):1.97e-2,
                 # str(pigment_dir+'ppg_shifted.csv'):1.97e-2,
                  #str(pigment_dir+'phenols_original_corrected.csv'):0.0,
                  #str(pigment_dir+'PhenolMAC_200nm.csv'):0.0,
                  str(pigment_dir+'phenols_invivo_corrected.csv'):2.53e-4,
                  #
                  }
for key,value in pigments_data.items():
            abs_pigm=(np.array(pd.read_csv(key))).flatten() # m2/mg_pigm 
            abs_cff_pigments[str(key.split(pigment_dir,1)[1])[0:-4]]=abs_pigm 
abs_cff_pigments.columns=['ppg standard', 'ppg from raw extract',  'phenols Williamson 2020', 'phenols unpubl. data C. Williamson']#, 'ppg standard shifted','phenols corrected Williamson 2020']

ax=abs_cff_pigments.replace(0, np.nan).plot(use_index=True,xlim=(270, 750), title="Comparison phenol/ppg absorption coefficient")
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("absorption coefficient (m$^2$ mg$^{-1}$)")
#ax.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/comp_ppg_cff_1.eps', format='eps')


#%%
pigments_data_glacier_algae_MIT5b_extract = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.61e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):4.66e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.04e-4,
                  str(pigment_dir+'pheophytin.csv'):3.61e-5,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.97e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):2.53e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.78e-4,
                  }

pigments_data_glacier_algae_MIT5b_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.61e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):4.66e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.04e-4,
                  str(pigment_dir+'pheophytin.csv'):3.61e-5,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_shifted.csv'):1.97e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):2.53e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.78e-4,
                  }

pigments_data_glacier_algae_MIT3a_extract = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.06e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):7.76e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.29e-5,
                  str(pigment_dir+'pheophytin.csv'):5.59e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):4.21e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.56e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT3a_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.06e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):7.76e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.29e-5,
                  str(pigment_dir+'pheophytin.csv'):5.59e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_shifted.csv'):4.21e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.56e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }


wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_extract_MIT5b, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_MIT5b_extract,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_reconstruction_MIT5b, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_MIT5b_shifted,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
    
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_extract_MIT3a, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_MIT3a_extract,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)
wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, MAC, MAC_rescaled_BioSNICAR_reconstruction_MIT3a, data, abs_coeff \
= bioptical_calculations(MAC_calculated, MAC_path, biovolume, density, xw, cell_volume, wvl, \
                       packaging_correction,pigment_dir, pigments_data_glacier_algae_MIT3a_shifted,\
                       n_algae, k_water, smooth, window_size,\
                       poly_order, smoothStart, smoothStop, plot_n_k_MAC_cell,\
                       plot_pigment_MACs, savefiles_n_k_MAC_cell, \
                       savefilename_n_k_MAC_cell, saveplots_n_k_MAC,\
                      savepath_n_k_MAC)

MAC_plot=pd.DataFrame(index=wvl[::10])
test=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/BioOptical_Model/algae_geom_4_40.nc')
MAC_plot['Cook 2020']=np.concatenate([np.zeros(10), np.array(test['ext_cff_mss']/(1*10**12))])
#MAC_plot['Raw extract Halbach 2021, MIT5b']=MAC_rescaled_BioSNICAR_extract_MIT5b
#MAC_plot['Reconstruction Halbach 2021, MIT5b']=MAC_rescaled_BioSNICAR_reconstruction_MIT5b
#MAC_Canada=xr.open_dataset('/Users/au660413/Desktop/GitHub/BioSNICAR_GO_PY/Data/Mie_files/480band/lap/Glacier_Algae_480.nc')
#MAC_plot['Chris_MAC_Canada']=(MAC_Canada['ext_cff_mss'])*1*10**(-12)
MAC_plot['Raw extract Halbach 2021, MIT3a']=MAC_rescaled_BioSNICAR_extract_MIT3a
MAC_plot['Reconstruction Halbach 2021, MIT3a']=MAC_rescaled_BioSNICAR_reconstruction_MIT3a
ax = MAC_plot.plot(use_index=True, xlim=(0.3,0.8), title='Glacier algae cellular MAC - Sermilik 2019')
ax.set_xlabel("wvl (µm)")
ax.set_ylabel("MAC (m$^2$ cell$^{-1}$)")
ax.figure.savefig('/Users/au660413/Documents/Aarhus_PhD/Data/Halbach_2021/raw_extract_MIT3b_baselined.eps', format='eps')

#%% 
#########################
## PIGMENT DATA 
#########################
##### COMPARISON PHENOL DATA FROM CHRIS 
pigments_data_allpig = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):4.7e-4,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):1.5e-3,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):1.6e-4,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):2.8e-3,
                  str(pigment_dir+'PhenolMAC_200nm.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):8e-5,
                  str(pigment_dir+'violaxanthin.csv'):5.7e-4,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):5.3e-4,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_phenols_original = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'phenols_original.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_phenols_corrected = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'phenols_original_corrected.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_phenols_invivocorrected = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'phenols_invivo_corrected.csv'):1e-2,#4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_ppg = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'ppg.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_ppg_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'ppg_shifted.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }
pigments_data_ppg_extract = {str(pigment_dir+'alloxanthin.csv'):0.0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#max=0,#mean=1.05e-11,min=0
                  str(pigment_dir+'chl-a.csv'):3.96e-3,#max=2.92e-9,#mean=1.38e-9,min=1.07e-10
                  str(pigment_dir+'chl-b.csv'):7e-4,#max=5.32e-10,#mean=3.24e-10,min=1.09e-9
                  str(pigment_dir+'lutein.csv'):0,#max=9.22e-10,#mean=2.6e-10,min=0
                  str(pigment_dir+'neoxanthin.csv'):0,#max=1.09e-10,#mean=3.21e-10,min=0
                  str(pigment_dir+'pheophytin.csv'):0.0,#max=0,#mean=1.22e-10,min=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,#1.3e-3, #echinone + cantaxanthin
                  str(pigment_dir+'Photos_carotenoids.csv'):6e-3,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):4.3e-2,#max=5.65e-8,#mean=2.45e-8, min=2.61e-9
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):0,#max=2.07e-10,#mean=3.02e-10,#min=2.52e-11
                  str(pigment_dir+'zeaxanthin.csv'):0,#max=1.90e-10,#mean=4.79e-11,min=0
                  }

# STATIONS FROM LAURA - ONLY ALGAE

pigments_data_glacier_algae_BR1_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.08e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):2.97e-5,
                  str(pigment_dir+'pheophytin.csv'):3.47e-3,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):3.68e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):7.06e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_BR2_shifted = {str(pigment_dir+'alloxanthin.csv'):3.68e-5,
                  str(pigment_dir+'antheraxanthin.csv'):3.68e-5,
                  str(pigment_dir+'chl-a.csv'):4.03e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1.02e-3,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.23e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):5.13e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_BR3_shifted = {str(pigment_dir+'alloxanthin.csv'):0.685e-4,
                  str(pigment_dir+'antheraxanthin.csv'):0.685e-4,
                  str(pigment_dir+'chl-a.csv'):4.47e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):1.90e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.87e-3,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.43e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):1.07e-3,
                  str(pigment_dir+'zeaxanthin.csv'):5.96e-4,
                  }

pigments_data_glacier_algae_HE2a_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.39e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0,
                  str(pigment_dir+'neoxanthin.csv'):1.30e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.66e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.78e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_HE2b_shifted = {str(pigment_dir+'alloxanthin.csv'):4.07e-5,
                  str(pigment_dir+'antheraxanthin.csv'):4.07e-5,
                  str(pigment_dir+'chl-a.csv'):8.93e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):1.90e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.02e-4,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.72e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.53e-4,
                  str(pigment_dir+'zeaxanthin.csv'):9.52e-6,
                  }

pigments_data_glacier_algae_HE2c_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.20e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):9.28e-6,
                  str(pigment_dir+'neoxanthin.csv'):3.86e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.46e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.21e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT6_shifted = {str(pigment_dir+'alloxanthin.csv'):5e-7,
                  str(pigment_dir+'antheraxanthin.csv'):5e-7,
                  str(pigment_dir+'chl-a.csv'):1.54e-4,
                  str(pigment_dir+'chl-b.csv'):4.90e-5,
                  str(pigment_dir+'lutein.csv'):3.15e-5,
                  str(pigment_dir+'neoxanthin.csv'):4.46e-6,
                  str(pigment_dir+'pheophytin.csv'):1.39e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):7.97e-3,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0,
                  str(pigment_dir+'violaxanthin.csv'):1.90e-5,
                  str(pigment_dir+'zeaxanthin.csv'):7.95e-6,
                  }
#######
pigments_data_glacier_algae_MIT3a_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.06e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):7.76e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.29e-5,
                  str(pigment_dir+'pheophytin.csv'):5.59e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):4.21e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):1.56e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_MIT3b_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.92e-3,
                  str(pigment_dir+'chl-b.csv'):5.32e-4,
                  str(pigment_dir+'lutein.csv'):9.22e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.09e-4,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.65e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):2.07e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.97e-4,
                  }
#####
pigments_data_glacier_algae_MIT3c_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.69e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):8.45e-5,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):1.71e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.19e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):4.94e-5,
                  str(pigment_dir+'zeaxanthin.csv'):8.43e-6,
                  }

pigments_data_glacier_algae_MIT4a_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.05e-3,
                  str(pigment_dir+'chl-b.csv'):7.82e-4,
                  str(pigment_dir+'lutein.csv'):1.34e-3,
                  str(pigment_dir+'neoxanthin.csv'):9.30e-4,
                  str(pigment_dir+'pheophytin.csv'):2.24e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.46e-3,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):7.90e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }


pigments_data_glacier_algae_MIT4b_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.25e-3,
                  str(pigment_dir+'chl-b.csv'):7.04e-4,
                  str(pigment_dir+'lutein.csv'):6.66e-4,
                  str(pigment_dir+'neoxanthin.csv'):8.39e-4,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.37e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):9.31e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }


pigments_data_glacier_algae_MIT4c_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.48e-4,
                  str(pigment_dir+'chl-b.csv'):3.04e-4,
                  str(pigment_dir+'lutein.csv'):4.32e-4,
                  str(pigment_dir+'neoxanthin.csv'):4.22e-5,
                  str(pigment_dir+'pheophytin.csv'):1.99e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.70e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):3.22e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.41e-4,
                  }

pigments_data_glacier_algae_MIT5a_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.07e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.61e-3,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):2.52e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_MIT5b_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.61e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):4.66e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.04e-4,
                  str(pigment_dir+'pheophytin.csv'):3.61e-5,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.97e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):2.53e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.78e-4,
                  }
pigments_data_glacier_algae_MIT5c_shifted = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):9.47e-4,
                  str(pigment_dir+'chl-b.csv'):7.69e-4,
                  str(pigment_dir+'lutein.csv'):4.71e-4,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):1.15e-5,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):3.26e-2,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin.csv'):0.0,
                  str(pigment_dir+'violaxanthin.csv'):3.20e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.68e-4,
                  }


#### SNOW ALGAE ONLY
pigments_data_snow_algae_MIT2c = {str(pigment_dir+'alloxanthin.csv'):0,#min=0,#mean=0,max=0
                  str(pigment_dir+'antheraxanthin.csv'):0,#min=0,#mean=0,max=0
                  str(pigment_dir+'chl-a.csv'):9.99e-3,#min=2.92e-9,#mean=8.86e-9,max=2.50e-8
                  str(pigment_dir+'chl-b.csv'):0,#min=0,#mean=0,max=0
                  str(pigment_dir+'lutein.csv'):3.12e-3,#min=6.31e-10,#mean=7.64e-10,max=4.09e-10
                  str(pigment_dir+'neoxanthin.csv'):1.25e-3,#min=0,#mean=1.65e-10,max=0
                  str(pigment_dir+'pheophytin.csv'):6.99e-3,#min=0,#mean=2.27e-9,max=0
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  #str(pigment_dir+'trans_astaxanthin.csv'):4.49e-4,#min=1.02e-10,#mean=2.67e-10,max=8.4e-10
                  #str(pigment_dir+'cis_astaxanthin_monoester.csv'): 1.47e-1,#min=0.0414e-6,#mean=1.81e-7,max=5e-7
                  #str(pigment_dir+'cis_astaxanthin_diester.csv'): 3.6e-1,#min=0.0414e-6,#mean=1.81e-7,max=5e-7
                  #str(pigment_dir+'trans_astaxanthin_monoester.csv'): 2.83e-2,#min=0.0414e-6,#mean=1.81e-7,max=5e-7
                  #str(pigment_dir+'trans_astaxanthin_diester.csv'): 6.29e-2,#min=0.0414e-6,#mean=1.81e-7,max=5e-7
                  str(pigment_dir+'total_astaxanthin.csv'): 5.98e-1,#min=0.0414e-6,#mean=1.81e-7,max=5e-7
                  str(pigment_dir+'violaxanthin.csv'):1.01e-11,#min=0,#mean=1.01e-11, max=0
                  str(pigment_dir+'zeaxanthin.csv'):8.08e-4,#min=0,#mean=1.52e-10,max=0
                 }

#%% STATIONS FROM LAURA - COMMUNITY

pigments_data_glacier_algae_BR1_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.58e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1.06e-4,
                  str(pigment_dir+'pheophytin.csv'):3.70e-3,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):3.68e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.301,
                  str(pigment_dir+'violaxanthin.csv'):6.61e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_BR2_community = {str(pigment_dir+'alloxanthin.csv'):3.68e-5,
                  str(pigment_dir+'antheraxanthin.csv'):3.68e-5,
                  str(pigment_dir+'chl-a.csv'):4.64e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1.04e-3,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.23e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.438,
                  str(pigment_dir+'violaxanthin.csv'):5.04e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_BR3_community = {str(pigment_dir+'alloxanthin.csv'):0.685e-4,
                  str(pigment_dir+'antheraxanthin.csv'):0.685e-4,
                  str(pigment_dir+'chl-a.csv'):5.05e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1.80e-3,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.43e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.487,
                  str(pigment_dir+'violaxanthin.csv'):9.59e-4,
                  str(pigment_dir+'zeaxanthin.csv'):6.18e-4,
                  }

pigments_data_glacier_algae_HE1a_community = {str(pigment_dir+'alloxanthin.csv'):4.67e-6,
                  str(pigment_dir+'antheraxanthin.csv'):4.67e-6,
                  str(pigment_dir+'chl-a.csv'):7.56e-5,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):2.66e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.95e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.001,
                  str(pigment_dir+'violaxanthin.csv'):3.76e-5,
                  str(pigment_dir+'zeaxanthin.csv'):5.57e-6,
                  }

pigments_data_glacier_algae_HE1b_community = {str(pigment_dir+'alloxanthin.csv'):1.86e-6,
                  str(pigment_dir+'antheraxanthin.csv'):1.86e-6,
                  str(pigment_dir+'chl-a.csv'):0.0,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1.15e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.002,
                  str(pigment_dir+'violaxanthin.csv'):1.75e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_HE1c_community = {str(pigment_dir+'alloxanthin.csv'):2.715e-5,
                  str(pigment_dir+'antheraxanthin.csv'):2.715e-5,
                  str(pigment_dir+'chl-a.csv'):2.78e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):4.10e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.009,
                  str(pigment_dir+'violaxanthin.csv'):6.95e-5,
                  str(pigment_dir+'zeaxanthin.csv'):6.66e-6,
                  }

pigments_data_glacier_algae_HE2a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.39e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0,
                  str(pigment_dir+'neoxanthin.csv'):1.30e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.66e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.004,
                  str(pigment_dir+'violaxanthin.csv'):1.78e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_HE2b_community = {str(pigment_dir+'alloxanthin.csv'):4.07e-5,
                  str(pigment_dir+'antheraxanthin.csv'):4.07e-5,
                  str(pigment_dir+'chl-a.csv'):8.93e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):1.90e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.02e-4,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.72e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.01,
                  str(pigment_dir+'violaxanthin.csv'):1.53e-4,
                  str(pigment_dir+'zeaxanthin.csv'):9.52e-6,
                  }

pigments_data_glacier_algae_HE2c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.20e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):9.28e-6,
                  str(pigment_dir+'neoxanthin.csv'):3.86e-5,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.46e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.03,
                  str(pigment_dir+'violaxanthin.csv'):1.21e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT6_community = {str(pigment_dir+'alloxanthin.csv'):5.02e-6,
                  str(pigment_dir+'antheraxanthin.csv'):5.02e-6,
                  str(pigment_dir+'chl-a.csv'):1.54e-4,
                  str(pigment_dir+'chl-b.csv'):4.90e-5,
                  str(pigment_dir+'lutein.csv'):3.15e-5,
                  str(pigment_dir+'neoxanthin.csv'):4.46e-6,
                  str(pigment_dir+'pheophytin.csv'):1.39e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):7.97e-3,
                  str(pigment_dir+'total_astaxanthin.csv'):0.019,
                  str(pigment_dir+'violaxanthin.csv'):1.90e-5,
                  str(pigment_dir+'zeaxanthin.csv'):7.95e-6,
                  }

pigments_data_glacier_algae_MIT1a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.32e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):9.52e-5,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.109,
                  str(pigment_dir+'violaxanthin.csv'):0.0,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT1b_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):4.60e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):4.72e-4,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.141,
                  str(pigment_dir+'violaxanthin.csv'):0.0,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT1c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.92e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):6.31e-4,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):0.0,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'trans_astaxanthin_ester.csv'):0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.076,
                  str(pigment_dir+'violaxanthin.csv'):0.0,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT2a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):3.10e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):0.0,
                  str(pigment_dir+'neoxanthin.csv'):1e-4,
                  str(pigment_dir+'pheophytin.csv'):3.43e-3,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.218,
                  str(pigment_dir+'violaxanthin.csv'):1.94e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT2b_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.61e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):1.12e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.17e-4,
                  str(pigment_dir+'pheophytin.csv'):2.69e-3,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.177,
                  str(pigment_dir+'violaxanthin.csv'):2.93e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }

pigments_data_glacier_algae_MIT2c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):9.99e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):3.12e-3,
                  str(pigment_dir+'neoxanthin.csv'):1.25e-3,
                  str(pigment_dir+'pheophytin.csv'):6.99e-3,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):0.0,
                  str(pigment_dir+'total_astaxanthin.csv'):0.527,
                  str(pigment_dir+'violaxanthin.csv'):0.0,
                  str(pigment_dir+'zeaxanthin.csv'):8.08e-4,
                  }


pigments_data_glacier_algae_MIT3a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.79e-3,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):9.66e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.50e-4,
                  str(pigment_dir+'pheophytin.csv'):1.08e-3,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):4.21e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.07,
                  str(pigment_dir+'violaxanthin.csv'):1.43e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_MIT3b_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):3.77e-3,
                  str(pigment_dir+'chl-b.csv'):4.68e-4,
                  str(pigment_dir+'lutein.csv'):1.19e-3,
                  str(pigment_dir+'neoxanthin.csv'):2.46e-4,
                  str(pigment_dir+'pheophytin.csv'):5.36e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.65e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.127,
                  str(pigment_dir+'violaxanthin.csv'):1.82e-4,
                  str(pigment_dir+'zeaxanthin.csv'):2.70e-4,
                  }
#####
pigments_data_glacier_algae_MIT3c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):3.80e-4,
                  str(pigment_dir+'chl-b.csv'):0.0,
                  str(pigment_dir+'lutein.csv'):1.50e-4,
                  str(pigment_dir+'neoxanthin.csv'):1.90e-5,
                  str(pigment_dir+'pheophytin.csv'):3.18e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.19e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.017,
                  str(pigment_dir+'violaxanthin.csv'):4.84e-5,
                  str(pigment_dir+'zeaxanthin.csv'):2.57e-5,
                  }

pigments_data_glacier_algae_MIT4a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.64e-3,
                  str(pigment_dir+'chl-b.csv'):7.25e-4,
                  str(pigment_dir+'lutein.csv'):1.47e-3,
                  str(pigment_dir+'neoxanthin.csv'):9.53e-4,
                  str(pigment_dir+'pheophytin.csv'):7.22e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):5.46e-3,
                  str(pigment_dir+'total_astaxanthin.csv'):0.149,
                  str(pigment_dir+'violaxanthin.csv'):7.32e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }


pigments_data_glacier_algae_MIT4b_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.34e-3,
                  str(pigment_dir+'chl-b.csv'):6.17e-4,
                  str(pigment_dir+'lutein.csv'):9.72e-4,
                  str(pigment_dir+'neoxanthin.csv'):8.90e-4,
                  str(pigment_dir+'pheophytin.csv'):4.93e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.37e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.084,
                  str(pigment_dir+'violaxanthin.csv'):8.16e-4,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }


pigments_data_glacier_algae_MIT4c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):5.76e-4,
                  str(pigment_dir+'chl-b.csv'):2.94e-4,
                  str(pigment_dir+'lutein.csv'):5.22e-4,
                  str(pigment_dir+'neoxanthin.csv'):8.28e-5,
                  str(pigment_dir+'pheophytin.csv'):4.28e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.70e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.034,
                  str(pigment_dir+'violaxanthin.csv'):3.12e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.63e-4,
                  }

pigments_data_glacier_algae_MIT5a_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.14e-3,
                  str(pigment_dir+'chl-b.csv'):9.79e-4,
                  str(pigment_dir+'lutein.csv'):3.22e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.22e-5,
                  str(pigment_dir+'pheophytin.csv'):3.94e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):2.61e-3,
                  str(pigment_dir+'total_astaxanthin.csv'):0.03,
                  str(pigment_dir+'violaxanthin.csv'):2.25e-5,
                  str(pigment_dir+'zeaxanthin.csv'):0.0,
                  }
pigments_data_glacier_algae_MIT5b_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):2.08e-3,
                  str(pigment_dir+'chl-b.csv'):1.21e-3,
                  str(pigment_dir+'lutein.csv'):6.15e-4,
                  str(pigment_dir+'neoxanthin.csv'):5.46e-4,
                  str(pigment_dir+'pheophytin.csv'):4.25e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):1.97e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.102,
                  str(pigment_dir+'violaxanthin.csv'):2.39e-4,
                  str(pigment_dir+'zeaxanthin.csv'):2.14e-4,
                  }
pigments_data_glacier_algae_MIT5c_community = {str(pigment_dir+'alloxanthin.csv'):0.0,
                  str(pigment_dir+'antheraxanthin.csv'):0.0,
                  str(pigment_dir+'chl-a.csv'):1.16e-3,
                  str(pigment_dir+'chl-b.csv'):7.51e-4,
                  str(pigment_dir+'lutein.csv'):5.33e-4,
                  str(pigment_dir+'neoxanthin.csv'):0.0,
                  str(pigment_dir+'pheophytin.csv'):2.76e-4,
                  str(pigment_dir+'Photop_carotenoids.csv'):0.0,
                  str(pigment_dir+'Photos_carotenoids.csv'):0.0,
                  str(pigment_dir+'ppg_raw_extract_mean.csv'):3.26e-2,
                  str(pigment_dir+'total_astaxanthin.csv'):0.048,
                  str(pigment_dir+'violaxanthin.csv'):3.12e-4,
                  str(pigment_dir+'zeaxanthin.csv'):1.83e-4,
                  }


pigments_data_HE1=clt.Counter({})
pigments_data_HE2=clt.Counter({})
pigments_data_MIT1=clt.Counter({})
pigments_data_MIT2=clt.Counter({})
pigments_data_MIT3=clt.Counter({})
pigments_data_MIT4=clt.Counter({})
pigments_data_MIT5=clt.Counter({})
pigments_data_BR=clt.Counter({})
nb_stations=3
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_HE1a_community, pigments_data_glacier_algae_HE1b_community, pigments_data_glacier_algae_HE1c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_HE1[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_HE2a_community, pigments_data_glacier_algae_HE2b_community, pigments_data_glacier_algae_HE2c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_HE2[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_MIT1a_community, pigments_data_glacier_algae_MIT1b_community, pigments_data_glacier_algae_MIT1c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_MIT1[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_MIT2a_community, pigments_data_glacier_algae_MIT2b_community, pigments_data_glacier_algae_MIT2c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_MIT2[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_MIT3a_community, pigments_data_glacier_algae_MIT3b_community, pigments_data_glacier_algae_MIT3c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_MIT3[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_MIT4a_community, pigments_data_glacier_algae_MIT4b_community, pigments_data_glacier_algae_MIT4c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_MIT4[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_MIT5a_community, pigments_data_glacier_algae_MIT5b_community, pigments_data_glacier_algae_MIT5c_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_MIT5[key]=total[key]/nb_stations
total=clt.Counter({})
for dictionary in [pigments_data_glacier_algae_BR1_community, pigments_data_glacier_algae_BR2_community, pigments_data_glacier_algae_BR3_community]:
    total=total+clt.Counter(dictionary)
for key in total:
    pigments_data_BR[key]=total[key]/nb_stations
    
    
    
    
    