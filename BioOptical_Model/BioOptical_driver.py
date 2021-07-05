"""
Created on Fri Jul 20 15:24:05 2018

@author: joe

This code supersedes Algae_GO.py, Algae_MIE.py and BioOptical_Model.py in the previous BioSNICAR_GO_PY commit, and
adds a driver script. This format is much more understandable and easier to handle.

The file BioOptical_Funcs now contains all the necessary functions for determining the optical properties
of algal cells from any of: a) measured cell MACs, b) pigment masses (mg), c) pigment mass fractions (%).
The result is the refractive index and mass absorption coefficient of cylindrical algal cells in a NetCDF
folder that can be accessed by the radiative transfer model BioSNICAR_GO to determine the algal impact
on snow or ice energy transfers. The pipeline flows as follows:

bio_optical() --> preprocess_RI --> calc_optical_params --> net_cdf_updater

In bio-optical() the MAC and/or imaginary refractive index of the cell is calculated. In preprocess_RI() the
relevant optical data is preprocessed to ensure equal length vectors with the correct resolution and wavelength range.
calc_optical_params uses the optical data from the first function along with the cell dimensions to calculate the 
single scattering optical properties of the cell, which are then passed to net_cdf_updater() for formatting and 
saving to the NetCDF library accessed by the radiative transfer model.

1) bio_optical()

There are several configuration options.

The user must first select whether the mass absorption coefficient of the algal cell is to be loaded in directly from a
single file (which would be the case if the user has made empirical measurements of cell MACs) or whether, instead the
cell MAC should be calculated from the amount of each pigment present in the cell. This is controlled by toggling
calc_MAC and load_MAC to True/False.

In the former case, the remaining parameters (apart from paths and plotting option) can be left as their defaults 
as they are not used in any further calculations. In the latter case, the user must define the units for the 
pigments. They can be provided as absolute masses of each pigment, in which case they should be given in mg, or 
they can be provided as mass fractions, in which case they should be values between 0-1 representing the fraction of
the total cellular dry mass comprised by that pigment. The user should toggle the parameters pig_mass and pig_frac 
to True/False accordingly. 

There is also an option to apply a correction for packaging effects in the cell. This is usually not applicable if
load_MAC = True, since empirical measurements of complete cells will naturally already accout for pigment packaging; 
however, it is recommended to apply this correction if building from individual pigment data. Setting 
apply_packaging_correction = True loads in a csv file containing the multiplier to be applied to each pigment's 
MAC at each unit wavelength.

The function will then return the cell MAC (m2/kg) and refractive index (dimensionless). The refractive index is
predicted using a mixing model adapted from either Pottier et al (2005) or Cook et al (2017). The user can select which
method is used by toggling Pottier and Cook to True/False. The benefit of Cook et al's (2017) method is that cells
have refractive indices identical to water at wavelengths > 750 nm rather than being completely non-absorbing, and there 
is no chance of the refractive index going infinite.

INPUT PARAMETERS:

load_MAC: if True, MAC loaded from external file
apply_packaging_correction: correct the MAC for packaging effects using a wavelength dependent correction factor (CW)
calc_MAC: If true, MAC calculated theoretically
calc_k: if True, script calculates refractive index
pig_mass: if True, pigment data should be provided as absolute mass of pigment per cell (mg)
pig_frac: if True, pigment data should be provided as mass fraction of pigment per cell (% total dry weight)
Pottier: If true, use Pottier et al's (2005) equation to predict the imaginary refractive index
Cook: If true, use Cook et al's (2019) updated equation to predict the imaginary refractive index
cell_dm_weight: provide dry weight of cell in ng
chla,chlb,ppro,psyn,purp: either mass (mg) or mass fraction of each pigment in cell
Xw: water fraction, set to 0.8 as default
density: density of dry cellular material in kg/m3, set to 1400 as default
nm: real part of refractive index, set to 1.4 as default
savefiles: if True, the data will be saved as a csv file to path specified in savefilename
savefilename: path to save files
plot_title: title for plots
plot_optical_props: if true, k and MAC plotted in ipython console
plot_pigment_MACs: if true, the mass absorption coefficient of each pigment is plotted to ipython console

OUTPUTS:

k_list: imaginary refractive index per unit wavelength from 300-5000 nm in at 10nm resolution
MAC: mass absorption coefficient per unit wavelength from 300-5000 nm at 10nm resolution
data: pandas dataframe containing nm, MAC, k per unit wavelength from 300 - 5000 nm at 10nm resolution

2) preprocess_RI()

This is a preprocessing function that ensures the wavelengths and real/imaginary parts of the refractive 
index for ice is provided in the correct waveband and correct spectral resolution to interface with the
BioSNICAR_GO model. 

3) calc_optical_params()

There are two flavours of this function - MIE and GO depending upon whether the algae are small and circular
in which cae MIE is appropriate or long and cylindrical in which case GO is appropriate.

This function calculates the optical properties (single scattering albedo, assymetry parameter, mass 
absorption coefficient and extinction, scattering and absorption cross sections) for algal cells shaped
as arbitrarily large cylinders. In GO mode, the main calculations are based upon the equations of 
Diedenhoven et al (2014) who provided a python script as supplementary material for their paper:

"A flexible parameterization for shortwave optical properties of ice crystals" by
Bastiaan van Diedenhoven; Andrew S. Ackerman; Brian Cairns; Ann M. Fridlind
accepted for publication in J. Atmos. Sci. (2013)

The original code can be downloaded from:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

in MIE mode, the optical properties are calculated using Mie scattering usign Scott prahl's miepython
package https://github.com/scottprahl/miepython. MiePython is used to calculate the single scattering 
optical properties from the refractive index and cell radius

calc_optical_params() takes several inputs: reals, imags and wavelengths as returned by preprocess_RI()
and user-defined values for cell radius and length. 

The final function, net_cdf_updater() is used to dump the optical parameters and
metadata into a netcdf file and save it into the working directory to be used as
a lookup library for the two-stream radiative transfer model BioSNICAR_GO.

NOTE: The extinction coefficient in the GO implementation == 2 for all size parameters
as assumed in the conventional geometric optics approximation.

"""


##################
# USER INPUTS
##################

from BioOptical_Funcs import bio_optical, preprocess_RI, calc_optical_params_GO, calc_optical_params_MIE, net_cdf_updater
import numpy as np
import pandas as pd

basePath = '/home/joe/Code/BioSNICAR_GO_PY/' # save path for full spectrum version
RealsPath = '/home/joe/Code/BioSNICAR_GO_PY/Data/Cell_InVivoPhenol_Real.csv'
ImagsPath = '/home/joe/Code/BioSNICAR_GO_PY/Data/Cell_InVivoPhenol_KK.csv'
MACPath = '/home/joe/Code/BioSNICAR_GO_PY/Data/Cell_InVivoPhenol_MAC.csv'
NetCDFpath = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/'
savepath=basePath




#####################
# FUNCTION CALLS
#####################

# 1) bio-optical()
k_list, real_list, MAC, data = bio_optical(basePath,load_MAC,\
        apply_packaging_correction, calc_MAC,calc_k,pig_mass,pig_frac,Pottier,Cook,\
        cell_dm_weight,chla,chlb,ppro,psyn, purp,Xw, cell_density, nm, smooth,\
        smoothStart,smoothStop,window_size,poly_order,savefiles,savepath,\
        savefilename = f"Cell_InVivoPhenol_{cell_radius}_{cell_length}_",\
        plot_optical_props,plot_pigment_MACs,saveplots)

# 2) preprocess_RI()
reals, imags, MAC, wavelengths = preprocess_RI(RealsPath,ImagsPath,MACPath)

# 3) calc_optical_params_GO() or calc_optical_params_MIE()
if GO:

    Assy_list, SSA_list, absXS_list, depth, r, \
        Chi_abs_list, Reff, X_list = calc_optical_params_GO(
            basePath, 
            cell_radius, 
            cell_length, 
            reals, 
            imags, 
            wavelengths, 
            plots=True, 
            report_dims=True)

elif MIE:
  
   Assy_list, SSA_list = calc_optical_params_MIE(
       basePath, 
       cell_radius, 
       cell_density, 
       reals, 
       imags, 
       wavelengths, 
       plots=True, 
       savefigs=False) 

# 4) net_cdf_updater()
net_cdf_updater(
    NetCDFpath, 
    Assy_list, 
    SSA_list, 
    MAC, 
    cell_length, 
    cell_radius, 
    cell_density)


# example loop for interating GO code over many cell dimensions:

# for r in np.arange(1,10,2):
#    for depth in np.arange(1,80,5):
#            reals, imags, MAC, wavelengths = preprocess_RI()
#            Assy_list,SSA_list,absXS_list,MAC_list,depth,r,Chi_abs_list,Reff,X_list = calc_optical_params(r,depth,reals,imags,wavelengths,plots=False,
#            report_dims = True)
#            net_cdf_updater(NetCDFpath,Assy_list,SSA_list,MAC,cell_length,cell_radius,cell_density)