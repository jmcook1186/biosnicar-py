"""
#########################
SNICAR_Inverter.py
Joseph Cook, October 2019
#########################

This script loads two LUTs from disk: clean_ice_LUT and dirty_ice_LUT.
When provided with a target spectrum (e.g. a field albedo spectrum or spectral albedo values from
drone or satellite image), a two-stage process quantifies the ice configuration and particle
type and concentration that best simulates it.

This is achieved by first isolating the NIR albedo and comparing this to all spectral albedos in the
clean ice LUT. The grain dimensions and ice density associated with the spectral albedo that matches 
the target spectrum most accurately (lowest wavelength-integrated absolute error) are retrieved.

These physical parameters are then used to isolate a subset of the dirty-ice LUT (i.e. only those
indices where grain dimensions and density == those retrieved in stage 1). The full target spectrum is
then compared to each spectral albedo in the sliced dirty ice LUT. The combination of LAP types and 
mass concentrations that best simulate the target spectrum, with the ice physical parameters held 
constant, are thereby retrieved.

The LUTs were generated previosuly using the script SNICAR_LUT_builder.py.


NOTE: Now vectorized (29th Oct 2019): 1000000 runs in 1min 59s (real)

NOTE: Currently relies upon full spectral resolution in target spectral albedo (i.e. 470 wavelengths)
To apply this to remote sensing data specific wavelengths/wavebands will be required

IMPORTANT: THIS VERSION TAKES INDIVIDUAL SPECTRA AND APPLIES A 2 STEP RETRIEVAL USING
THE NIR AND THEN THE VIS WAVEBANDS. IT IS NOT FEASIBLE TO SCALE THIS TO LARGE REMOTE
SENSING IMAGES BECAUSE A) IT RELIES ON FULL SPECTRAL RESOLUTION, B) IT IS COMPUTATIONALLY EXPENSIVE

A MODIFIED VERSION HAS BEEN INTEGRATED INTO COOK & TEDSTONE'S BIG ICE CLASSIFIER SCRIPT

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


##############################################################################
# DEFINE VARIABLES AND TARGET SPECTRA
#####################################


# define wavelength range for clean ice matching
# (should normally be restricted to NIR wavelengths with indexes 40-90)
wvl_start = 40
wvl_stop = 90

# load sample spectrum (eventually this will be reflectance values pixelwise in RS
# dataset)
spectrum = pd.read_csv('/home/joe/Desktop/TEST.csv',header=None)
spectrum = np.squeeze(np.array(spectrum))

# load LUTs from disk
cleanLUT = np.load('/home/joe/Desktop/clean_ice_LUT.npy')
dirtyLUT = np.load('/home/joe/Desktop/dirty_ice_LUT.npy')

# copy range of values used to generate LUT. THESE MUST MATCH THOSE USED TO BUILD LUT!!!
depths = [[5000,5000,5000,5000,5000],[6000,6000,6000,6000,6000],[7000,7000,7000,7000,7000],[8000,8000,8000,8000,8000],[9000,9000,9000,9000,9000],[10000,10000,10000,10000,10000],[12000,12000,12000,12000,12000],[15000,15000,15000,15000,15000]]
side_lengths = [[5000,5000,5000,5000,5000],[6000,6000,6000,6000,6000],[7000,7000,7000,7000,7000],[8000,8000,8000,8000,8000],[9000,9000,9000,9000,9000],[10000,10000,10000,10000,10000],[12000,12000,12000,12000,12000],[15000,15000,15000,15000,15000]]
densities = [[200,200,200,200,200],[300,300,300,300,300],[400,400,400,400,400],[500,500,500,500,500],[600,600,600,600,600],[700,700,700,700,700],[800,800,800,800,800],[900,900,900,900,900]]
dust = [[100000,0,0,0,0],[500000,0,0,0,0],[800000,0,0,0,0],[1000000,0,0,0,0]]
algae = [[100000,0,0,0,0],[500000,0,0,0,0],[800000,0,0,0,0],[1000000,0,0,0,0]]
wavelengths = np.arange(0.3,5,0.01)


##############################################################################
# DEFINE FUNCTIONS
###################

def retrieve_ice_params(cleanLUT, wvl_start, wvl_stop, spectrum, depths, side_lengths, densities, wavelengths, plot_comparison=False):

    # reduce spectrum to NIR wavelengths
    sample_spectrum = spectrum[wvl_start:wvl_stop]

    # reformat LUT
    LUTshape = cleanLUT.shape # get shape
    cleanLUT = cleanLUT.reshape(LUTshape[0]*LUTshape[1]*LUTshape[2],len(wavelengths)) # flatten to 2D array 
    cleanLUT = cleanLUT[:,wvl_start:wvl_stop] # reduce wavelengths to match sample spectrum
    nbr_spectra = len(cleanLUT[:,0]) # count columns

    # calculate error columnwise
    error = abs(np.subtract(cleanLUT,sample_spectrum)) # subtract sample from target with uFunc
    mean_error = np.mean(error,axis=1) # column means
    
    # find index of minimum error
    reshaped_LUT_idx = np.argmin(mean_error)

    # reverse engineer index of min error spectrum in original 6D array from reshaped 2D array
    nbr_depths = len(depths)
    nbr_sidelengths = len(side_lengths)
    nbr_densities = len(densities)
    original_LUT_idx = np.unravel_index(reshaped_LUT_idx, (nbr_depths,nbr_sidelengths,nbr_densities))

    # find variable values from idxs
    side_length = side_lengths[original_LUT_idx[0]]
    depth = depths[original_LUT_idx[1]]
    density = densities[original_LUT_idx[2]]

    if plot_comparison:
        idx1 = original_LUT_idx[0]
        idx2 = original_LUT_idx[1]
        idx3 = original_LUT_idx[2]

        plt.plot(wavelengths[wvl_start:wvl_stop],sample_spectrum), plt.plot(wavelengths[wvl_start:wvl_stop],cleanLUT[reshaped_LUT_idx,:])
        plt.ylabel("Albedo"),plt.xlabel("Wavelength (microns)")
        plt.title("Sample spectrum and closest match: NIR wavelengths only")

    return side_length, depth, density, original_LUT_idx



def retrieve_impurities(dirtyLUT, original_LUT_idx, spectrum, depths, side_lengths, densities, wavelengths, dust, algae, plot_comparison = True):

    dirtyLUT = dirtyLUT[original_LUT_idx[0],original_LUT_idx[1],original_LUT_idx[2],:,:,:]
    nbr_dust = len(dust)
    nbr_alg = len(algae)
    dirtyLUT = dirtyLUT.reshape(nbr_dust*nbr_alg,len(wavelengths))

    error = abs(np.subtract(dirtyLUT,spectrum))
    mean_error = np.mean(error,axis=1)

    reshaped_LUT_idx = np.argmin(mean_error)
    original_LUT_idx = np.unravel_index(reshaped_LUT_idx, (nbr_dust,nbr_alg))

    # find values from idxs
    dust_vals = dust[original_LUT_idx[0]]
    algae_vals = algae[original_LUT_idx[1]]

    if plot_comparison:
        idx1 = original_LUT_idx[0]
        idx2 = original_LUT_idx[1]

        plt.plot(wavelengths,spectrum), plt.plot(wavelengths,dirtyLUT[reshaped_LUT_idx,:])
        plt.ylabel("Albedo"),plt.xlabel("Wavelength (microns)")
        plt.title("Sample spectrum and closest match: NIR wavelengths only")

    return dust_vals, algae_vals


######################################################
# RUN FUNCTIONS
###############


for i in range(1000000):
    side_length, depth, density, original_LUT_idx = retrieve_ice_params(cleanLUT, wvl_start, wvl_stop, spectrum, depths, side_lengths, densities, wavelengths, plot_comparison=False)
    dust_vals, algae_vals = retrieve_impurities(dirtyLUT, original_LUT_idx, spectrum, depths, side_lengths, densities, wavelengths, dust, algae, plot_comparison = False)

print("\n\nSNICAR CONFIGURATION THAT BEST MATCHES SAMPLE SPECTRUM: \n")
print(f"Grain side Lengths (microns) = {side_length}")
print(f"Grain depths (microns) = = {depth}")
print(f"Layer densities (kg m-3) = {density}")
print(f"Dust in upper 1mm (ppb) = {dust_vals}")
print(f"Algae in upper 1mm (ppb) = {algae_vals}")
print()
print()