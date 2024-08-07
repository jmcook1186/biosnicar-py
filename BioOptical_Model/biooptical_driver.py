#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Author: Joseph Cook (original writing); Lou-Anne Chevrollier

Driver for the bio-optical model developed by Cook et al. 2017, 2020
and Chevrollier et al. 2022 to calculate algae optical properties 
and save them to a netcdf files directly usable in BioSNICAR. T
he model workflow consists in three functions detailed
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
Outputs of the different functions
##############################################################################
    1. numpy arrays for wvl, n, k and absorption cross section
        in the spectral range 200-5000µm w/
        both 1nm resolution and 10nm (=BioSNICAR) resolution (xxx_rescaled)
        for a given cell;
        dataframe combining variables at BioSNICAR (10nm) resolution;
        dataframe with individual pigment absorption coefficients
    2. numpy arrays with ssa and g for a single cell in the spectral
        range 200-5000µm at BioSNICAR (10nm) resolution
    3. netcdf file with extinction cross section, ssa and g that 
        can be used directly in BioSNICAR


Note:
The main calculations used for the GO mode are based upon the equations of
Diedenhoven et al (2014), who provided a python script as supplementary
material for their paper, available at:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

In Mie mode, the optical properties are calculated using Mie scattering
using Scott Prahl's miepython package https://github.com/scottprahl/miepython.
"""

import numpy as np

from biooptical_Funcs import (
    gaussian,
    bioptical_calculations,
    net_cdf_updater,
    ssp_calculations,
)

# %%

# --------------------------------------------------------------------------------------
# INPUTS TO FILL
# --------------------------------------------------------------------------------------

# Set base directory and constant variables
DIR_BASE = "./biosnicar-py"
WVL = np.arange(0.200, 5, 0.001)  # spectral range of interest in µm

# Chose ACS units and k correction
BIOVOLUME = False
BIOMASS = False
CELLULAR = False
BASELINE_CORRECTION_K_ALGAE = False

# Chose if ACS is calculated from pigment abs coeff or loaded
ACS_LOADED_INVIVO = False
ACS_LOADED_RECONSTRUCTED = False
ACS_CALCULATED = False

# if ACS is loaded:
ACS_FILE = ""


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


# Medium properties
N_MEDIUM = np.ones(480)
K_WATER = 0 * np.ones(np.size(WVL))  # only used if biovolume-prescribed ACS

# Algae properties
N_ALGAE = 1.38 * np.ones(np.size(WVL))
DENSITY_DRY = 600
DENSITY_WET = 1060

# Set up normal distribution around R to average Mie resonance features
R = 15
L = 20
mu = R
sigma = 0.1 * R
radii = np.linspace(R - 4 * sigma, R + 4 * sigma, 200)
cell_volumes = [4 / 3 * np.pi * r**3 for r in radii]
normal_distribution = gaussian(radii, R, sigma)


# Chose method for calculation of scattering optical properties
GO = False
MIE = True

# Optional smoothing filter for calculated k, ACS
SMOOTH = False
WINDOW_SIZE = 8
POLY_ORDER = 1
SMOOTH_START = 0
SMOOTH_STOP = 1000

# Directories and printing/saving options for calculated k, ACS
PLOT_N_K_ACS_CELL = False
PLOT_PIGMENT_ACSS = False
SAVEFILES_N_K_ACS_CELL = False
SAVEPLOTS_N_K_ACS = False
SAVEPATH_N_K_ACS_PLOTS = DIR_BASE
SAVEFILENAME_N_K_ACS_CELL = ""

# Directories and printing/saving options for scattering OPs
PLOTS_OPS = False
SAVEFIGS_OPS = False
REPORT_DIMS = False
SAVEPATH_OPS = DIR_BASE
FIGNAME_OPS = "figname"

# Saving OPs in netcdf
NETCDF_SAVE = False
SAVEPATH_NETCDF = "./biosnicar-py/Data/OP_data/480band/lap/"
FILENAME_NETCDF = ""
INFO = ""

# %%

asm_prm, ss_alb, acs, qext = [], [], [], []

for i, radius, cell_vol in zip(np.arange(0, len(radii), 1), radii, cell_volumes):

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
        BASELINE_CORRECTION_K_ALGAE,
        DENSITY_WET,
        DENSITY_DRY,
        DIR_BASE,
        cell_vol,
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

    #
    # --------------------------------------------------------------------------------------
    # CALCULATIONS OF SCATTERING PROPERTIES
    # --------------------------------------------------------------------------------------

    assym, ssa, q_ext = ssp_calculations(
        GO,
        MIE,
        SAVEPATH_OPS,
        radius,
        L,
        WVL,
        N_MEDIUM,
        n,
        k,
        PLOTS_OPS,
        SAVEFIGS_OPS,
        FIGNAME_OPS,
        REPORT_DIMS,
    )

    asm_prm.append(assym)
    ss_alb.append(ssa)
    qext.append(q_ext)
    acs.append(ACS)


asm_prm_weighed = np.sum(np.array(asm_prm).T * normal_distribution, axis=1) / sum(
    normal_distribution
)
ss_alb_weighed = np.sum(np.array(ss_alb).T * normal_distribution, axis=1) / sum(
    normal_distribution
)
qext_weighed = (
    np.sum(np.array(qext).T * normal_distribution, axis=1)
    / sum(normal_distribution)
    * np.pi
    * (R * 1e-6) ** 2
)
acs_weighed = np.sum(np.array(acs).T * normal_distribution, axis=1) / sum(
    normal_distribution
)

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
        asm_prm_weighed[5::10],
        ss_alb_weighed[5::10],
        qext_weighed[5::10],
        acs_weighed[5::10],
        L,
        R,
        DENSITY_WET,
        INFO,
    )
