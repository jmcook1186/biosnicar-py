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

#%%
from biooptical_Funcs import bioptical_calculations, net_cdf_updater, ssp_calculations

#%%

# --------------------------------------------------------------------------------------
# INPUTS TO FILL
# --------------------------------------------------------------------------------------

# Set base directory and constant variables
dir_base = "path to package" + "/BioSNICAR_GO_PY/"
wvl = np.arange(0.200, 4.999, 0.001)  # spectral range of interest in µm
k_water = np.loadtxt(dir_base + "Data/OP_data/k_ice_480.csv")

# Chose ACS units
biovolume = False
biomass = False
cellular = True

## Chose if ACS is calculated from pigment abs coeff or loaded
ACS_loaded_invivo = True
ACS_loaded_reconstructed = False
ACS_calculated = False
# if ACS is loaded:
ACS_file = "filename"
# if reconstructed ACS is directly loaded from pigment absorbance:
packaging_correction_SA = False
packaging_correction_GA = True
# if ACS is reconstructed from pigment profiles:
pigment_dir = dir_base + "Data/pigments/"
pigments_data = {
    str(pigment_dir + "alloxanthin.csv"): 0.0,
    str(pigment_dir + "antheraxanthin.csv"): 0,
    str(pigment_dir + "chl-a.csv"): 3.96e-3,
    str(pigment_dir + "chl-b.csv"): 7e-4,
    str(pigment_dir + "lutein.csv"): 0,
    str(pigment_dir + "neoxanthin.csv"): 0,
    str(pigment_dir + "pheophytin.csv"): 0.0,
    str(pigment_dir + "Photop_carotenoids.csv"): 0.0,
    str(pigment_dir + "Photos_carotenoids.csv"): 6e-3,
    str(pigment_dir + "ppg_shifted.csv"): 4.3e-2,
    str(pigment_dir + "trans_astaxanthin_ester.csv"): 0.0,
    str(pigment_dir + "trans_astaxanthin.csv"): 0,
    str(pigment_dir + "violaxanthin.csv"): 0,
    str(pigment_dir + "zeaxanthin.csv"): 0,
}

# Algae properties
n_algae = 1.4 * np.ones(np.size(wvl))
r = 5
L = 20
cell_volume = 1500
density_dry = 684
density_wet = 1160

# Chose method for calculation of scattering optical properties
GO = True
Mie = False

# Optional smoothing filter for calculated k, ACS
smooth = True
window_size = 25
poly_order = 3
smoothStart = 44
smoothStop = 100

# Directories and printing/saving options for calculated k, ACS
plot_n_k_ACS_cell = True
plot_pigment_ACSs = False
savefiles_n_k_ACS_cell = False
saveplots_n_k_ACS = False
savepath_n_k_ACS_plots = dir_base
savefilename_n_k_ACS_cell = ""

# Directories and printing/saving options for scattering OPs
plots_OPs = True
savefigs_OPs = False
report_dims = False
savepath_OPs = dir_base
figname_OPs = "figname"

# Saving OPs in netcdf
netcdf_save = False
savepath_netcdf = dir_base + "Data/OP_data/480band/lap/"
filename_netcdf = "filename"
info = ""


#%%
# --------------------------------------------------------------------------------------
# CALCULATIONS OF ABSORPTION PROPERTIES
# --------------------------------------------------------------------------------------

(
    wvl,
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
    ACS_calculated,
    ACS_file,
    ACS_loaded_invivo,
    ACS_loaded_reconstructed,
    biovolume,
    biomass,
    cellular,
    density_wet,
    density_dry,
    dir_base,
    cell_volume,
    wvl,
    packaging_correction_SA,
    packaging_correction_GA,
    pigment_dir,
    pigments_data,
    n_algae,
    k_water,
    smooth,
    window_size,
    poly_order,
    smoothStart,
    smoothStop,
    plot_n_k_ACS_cell,
    plot_pigment_ACSs,
    savefiles_n_k_ACS_cell,
    savefilename_n_k_ACS_cell,
    saveplots_n_k_ACS,
    savepath_n_k_ACS_plots,
)


#%%
# --------------------------------------------------------------------------------------
# CALCULATIONS OF SCATTERING PROPERTIES
# --------------------------------------------------------------------------------------

assym, ss_alb = ssp_calculations(
    GO,
    Mie,
    savepath_OPs,
    r,
    L,
    wvl_rescaled_BioSNICAR,
    n_rescaled_BioSNICAR,
    k_rescaled_BioSNICAR,
    plots_OPs,
    savefigs_OPs,
    figname_OPs,
    report_dims,
)

#%%
# --------------------------------------------------------------------------------------
# SAVING DATA IN NETCDF
# --------------------------------------------------------------------------------------

if netcdf_save:
    net_cdf_updater(
        GO,
        Mie,
        savepath_netcdf,
        filename_netcdf,
        wvl_rescaled_BioSNICAR,
        assym,
        ss_alb,
        ACS_rescaled_BioSNICAR,
        L,
        r,
        density_wet,
        info,
    )
