#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joe, lou

Driver for the bio-optical model developed by Cook et al. 2017, 2020 to
calculate algae optical properties and save them to a netcdf files directly
usable in BioSNICAR. The model workflow consists in three functions detailed
below in case the user needs to use the functions separately.
If the user wants to run the full model, the parameters to document
are summarized at the beginning of the code.

##############################################################################
Model workflow
##############################################################################

    1. bioptical_calculations() : Calculation/loading of absorption cross
        section (ACS) in m2/cell, m2/um3 or m2/ng and calculation of
        refractive index (n,k; unitless) in the spectral range of interest,
        both at 1nm and 10nm resolution (used in BioSNICAR)
    2. ssp_calculations() : Calculation of single scattering properties
        (g, ssa; unitless) using Mie or geometric optics theory in the
        spectral range of interest at the resolution of BioSNICAR (10nm)
    3. save_netcdf() : Storage of the data and associated metadata in a netcdf
        file directly usable in BioSNICAR

##############################################################################
Inputs of the different functions
##############################################################################

    1.  wvl = (numpy array, default = np.arange(0.200, 4.999, 0.001))
            wavelengths in spectral range of interest (in µm, 1nm step)
        cellular = (boolean) ACS units of m2/cell (recommended)
        biovolume = (boolean) ACS units of m2/um3
        biomass = (boolean) ACS units of m2/mg
        ACS_calculated = (boolan) True if the ACS is
                        calculated from intracellular pigment concentrations,
                        False if loaded
        pigment_dir = (string) used if ACS_calculated is True, directory to
                     folder containing pigment mass absorption coefficients
                     that must be csv file with size and resolution of wvl,
                     and units in m2/mg
        pigment_data = (string) dictionary with pigment file names and
                        associated intracellular concentrations
                        (ng/cell, ng/µm3 or ng/ng)
        ACS_loaded_reconstructed = (boolean) True if the
                            ACS is loaded as a reconstructed spectrum
                            from pigment absorbance (see methods in
                            Chevrollier et al. 2022)
        packaging_correction = (boolean - applied ONLY if
                            ACS_loaded_reconstructed is True) if True,
                            reconstructed ACS is corrected for pigment
                            packaging following Chevrollier et al. 2022
        ACS_loaded_invivo = (boolean) True if the ACS is loaded as in vivo
                            spectra of whole cells
        ACS_file = (string) directory to the ACS file if loaded
        density_dry =  (int - used if biomass = True) density of dry biomass
                    (kg/m3 - 625 and 684 for snow and glacier algae,
                    Chevrollier et al. 2022)
        density_wet =  (int - used if biomass = True) density of wet biomass
                    (kg/m3 - 1060 and 1160 for snow and glacier algae,
                    Chevrollier et al. 2022)
        cell_volume = (int - used if cellular = True) volume of the algae
                    cell (um3)
        n_algae = (numpy array) real part of cellular refractive index
                    in the spectral range of wvl (constant 1.38 by default,
                    Chevrollier et al. 2022)
        k_water = (numpy array) imaginary part of the refractive index of water
                    in the spectral range of wvl
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

     2. GO = (boolean) if True, uses geometric optics equations (Cook et
        al. 2020 adapted from Diedenhoven et al (2014)) to calculate single
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
        info = (string) containing any additional info for
                    metadata in netcdf (e.g. 'Glacier algae OPs derived
                    from GO calculations with empirical ACS')

##############################################################################
Outputs of the different functions
##############################################################################
    1. numpy arrays for wvl, n, k and CAC in the spectral range 200-5000µm w/
        both 1nm resolution and 10nm (=BioSNICAR) resolution (xxx_rescaled)
        for a given cell;
        dataframe combining variables at BioSNICAR (10nm) resolution;
        dataframe with individual pigment absorption coefficients
    2. numpy arrays with ssa and g for a single cell in the spectral
        range 200-5000µm at BioSNICAR (10nm) resolution
    3. netcdf file with CAC, ssa and g that can be used directly in BioSNICAR


Note:
The main calculations used for the GO mode are based upon the equations of
Diedenhoven et al (2014), who provided a python script as supplementary
material for their paper, available at:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

In Mie mode, the optical properties are calculated using Mie scattering
using Scott Prahl's miepython package https://github.com/scottprahl/miepython.
"""
import numpy as np

# %%
from biooptical_Funcs import bioptical_calculations, net_cdf_updater, ssp_calculations

# %%

# --------------------------------------------------------------------------------------
# INPUTS TO FILL
# --------------------------------------------------------------------------------------

# Set base directory and constant variables
DIR_BASE = "./BioSNICAR_GO_PY/"
WVL = np.arange(0.200, 4.999, 0.001)  # spectral range of interest in µm
K_WATER = 0 * np.ones(np.size(WVL)) #np.loadtxt(DIR_BASE + "Data/OP_data/k_ice_480.csv")

# Chose ACS units
BIOVOLUME = False
BIOMASS = True
CELLULAR = False

# Chose if ACS is calculated from pigment abs coeff or loaded
ACS_LOADED_INVIVO = True
ACS_LOADED_RECONSTRUCTED = False
ACS_CALCULATED = False
# if ACS is loaded:
ACS_FILE = "/Users/au660413/Documents/AU_PhD/ms2_dust_albedo/icesurf_MAC_filter_m2mg.csv"
# if reconstructed ACS is directly loaded from pigment absorbance:
PACKAGING_CORRECTION_SA = False
PACKAGING_CORRECTION_GA = False
# if ACS is reconstructed from pigment profiles:
PIGMENT_DIR = DIR_BASE + "Data/pigments/"
PIGMENTS_DATA = {
    str(PIGMENT_DIR + "alloxanthin.csv"): 0.0,
    str(PIGMENT_DIR + "antheraxanthin.csv"): 0,
    str(PIGMENT_DIR + "chl-a.csv"): 3.96e-3,
    str(PIGMENT_DIR + "chl-b.csv"): 7e-4,
    str(PIGMENT_DIR + "lutein.csv"): 0,
    str(PIGMENT_DIR + "neoxanthin.csv"): 0,
    str(PIGMENT_DIR + "pheophytin.csv"): 0.0,
    str(PIGMENT_DIR + "Photop_carotenoids.csv"): 0.0,
    str(PIGMENT_DIR + "Photos_carotenoids.csv"): 6e-3,
    str(PIGMENT_DIR + "ppg_shifted.csv"): 4.3e-2,
    str(PIGMENT_DIR + "trans_astaxanthin_ester.csv"): 0.0,
    str(PIGMENT_DIR + "trans_astaxanthin.csv"): 0,
    str(PIGMENT_DIR + "violaxanthin.csv"): 0,
    str(PIGMENT_DIR + "zeaxanthin.csv"): 0,
}

# Algae properties
N_ALGAE = 1.55 * np.ones(np.size(WVL))
R = [0.45]
L = 2.4
CELL_VOLUME = 4/3*np.pi*R[0]**3 #L*np.pi*R[0]**2 
DENSITY_DRY = 2500
DENSITY_WET = 0

# Chose method for calculation of scattering optical properties
GO = False
MIE = True

# Optional smoothing filter for calculated k, ACS
SMOOTH = True
WINDOW_SIZE = 25
POLY_ORDER = 3
SMOOTH_START = 44
SMOOTH_STOP = 100

# Directories and printing/saving options for calculated k, ACS
PLOT_N_K_ACS_CELL = True
PLOT_PIGMENT_ACSS = False
SAVEFILES_N_K_ACS_CELL = False
SAVEPLOTS_N_K_ACS = False
SAVEPATH_N_K_ACS_PLOTS = DIR_BASE
SAVEFILENAME_N_K_ACS_CELL = ""

# Directories and printing/saving options for scattering OPs
PLOTS_OPS = True
SAVEFIGS_OPS = False
REPORT_DIMS = False
SAVEPATH_OPS = DIR_BASE
FIGNAME_OPS = "figname"

# Saving OPs in netcdf
NETCDF_SAVE = True
SAVEPATH_NETCDF = "/Users/au660413/Desktop/github/biosnicar-py/Data/OP_data/480band/lap/"
FILENAME_NETCDF = "dust_QQT_Mie_1.6um"
INFO = ""


# %%
    # --------------------------------------------------------------------------------------
    # CALCULATIONS OF ABSORPTION PROPERTIES
    # --------------------------------------------------------------------------------------
    
(
    WVL,
    wvl_rescaled_BioSNICAR,
    k,
    k_rescaled_BioSNICAR,
    n,
    n_rescaled_BioSNICAR,
    ACS,
    ACS_rescaled_BioSNICAR,
    n_k_ACS_rescaled_BioSNICAR,
    abs_coeff_pigm_DataFrame,
) = bioptical_calculations(
    ACS_CALCULATED,
    ACS_FILE,
    ACS_LOADED_INVIVO,
    ACS_LOADED_RECONSTRUCTED,
    BIOVOLUME,
    BIOMASS,
    CELLULAR,
    DENSITY_WET,
    DENSITY_DRY,
    DIR_BASE,
    CELL_VOLUME,
    WVL,
    PACKAGING_CORRECTION_SA,
    PACKAGING_CORRECTION_GA,
    PIGMENT_DIR,
    PIGMENTS_DATA,
    N_ALGAE,
    K_WATER,
    SMOOTH,
    WINDOW_SIZE,
    POLY_ORDER,
    SMOOTH_START,
    SMOOTH_STOP,
    PLOT_N_K_ACS_CELL,
    PLOT_PIGMENT_ACSS,
    SAVEFILES_N_K_ACS_CELL,
    SAVEFILENAME_N_K_ACS_CELL,
    SAVEPLOTS_N_K_ACS,
    SAVEPATH_N_K_ACS_PLOTS,
)

    
# %%
# --------------------------------------------------------------------------------------
# CALCULATIONS OF SCATTERING PROPERTIES
# --------------------------------------------------------------------------------------
# RAD = [0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.75,1.25,1.75,2.25,3.75,7.5,12.5,17.5]
# weights_nb = [3.282275711,10.94091904,13.12910284,10.50328228,
#            8.533916849,5.251641138,4.376367615,3.938730853,
#            2.84463895,2.188183807,13.78555799,6.564551422,
#            4.157549234,1.969365427,4.814004376,2.84463895,
#            0.656455142,0.218818381]
# weights_vol = [0.001649385,0.01650214,0.033054136,0.037114072,
#                0.038908435,0.0293985,0.029116059,0.030438415,
#                0.025108443,0.021778676,0.218818381,0.176466436,
#                0.160789197,0.101863729,0.437636761,0.56892779,
#                0.273522976,0.218818381]
# assym_final = np.zeros(np.size(k_rescaled_BioSNICAR))
# ssa_final = np.zeros(np.size(k_rescaled_BioSNICAR))
RAD=[1.6]
for i,radius in enumerate(RAD):
    
    assym, ss_alb = ssp_calculations(
        GO,
        MIE,
        SAVEPATH_OPS,
        radius,
        L,
        wvl_rescaled_BioSNICAR,
        n_rescaled_BioSNICAR,
        k_rescaled_BioSNICAR,
        PLOTS_OPS,
        SAVEFIGS_OPS,
        FIGNAME_OPS,
        REPORT_DIMS,
    )
    #assym_final = assym_final + assym*weights_vol[i]/100
    #ssa_final = ssa_final + ss_alb*weights_vol[i]/100

    

# %%
# --------------------------------------------------------------------------------------
# SAVING DATA IN NETCDF
# --------------------------------------------------------------------------------------

if NETCDF_SAVE:
    net_cdf_updater(
        GO,
        MIE,
        SAVEPATH_NETCDF,
        FILENAME_NETCDF,
        wvl_rescaled_BioSNICAR,
        assym,
        ss_alb,
        ACS_rescaled_BioSNICAR,
        L,
        1,
        DENSITY_WET,
        INFO,
    )
