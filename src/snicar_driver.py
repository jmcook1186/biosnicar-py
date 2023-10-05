"""
BioSNICAR_GO DRIVER SCRIPT

This script is used to configure the 2-stream radiative transfer
model BioSNICAR_GO. Here variable values are defined, the model called
and the results plotted.

NB. Setting Mie = 1, GO = 0 and algal impurities = 0 is equivalent to
running the original SNICAR model of Flanner et al. (2007, 2009)

NB: if using only granular layers, recommend using the faster toon et al
tridiagonal matix solver (by setting toon = True), however this will not
include any specular reflection components. If solid ice layers are
included in the ice/snow column, the ADDING-DOUBLING solver must be used
(i.e. add_double = True).

glacier algae MAC files are provided in units of m2/cell, which requires
a unit conversion that is not applied to the other LAPS. The conversion
is selectivlely applid by indexing to the last element in the LAP
lists, meaning glacier algae must always be the final LAP, even if more
LAPs are added in future.

Author: Joseph Cook, June 2021
"""


import collections as c
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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


# --------------------------------------------------------------------------------------
# 2) Set working directory
# --------------------------------------------------------------------------------------

# set dir_base to the location of the BioSNICAR_GO_PY folder
Inputs.dir_base = str(Path(__file__).parent.parent) + "/"
savepath = Inputs.dir_base  # base path for saving figures
WRITE_CONFIG_TO_TEXTFILE = False  # toggle to save config to file
Inputs.verbosity = 0  # 1 to print real-time updates


# --------------------------------------------------------------------------------------
# 3) Choose plot/print options
# --------------------------------------------------------------------------------------

SHOW_FIGS = True  # toggle to display spectral albedo figure
SAVE_FIGS = False  # toggle to save spectral albedo figure to file
PRINT_BBA = True  # toggle to print broadband albedo to terminal
PRINT_BAND_RATIOS = True  # toggle to print various band ratios to terminal
SMOOTH = False  # apply optional smoothing function (Savitzky-Golay filter)
WINDOW_SIZE = 9  # if applying smoothing filter, define window size
POLY_ORDER = 3  # if applying smooting filter, define order of polynomial

# --------------------------------------------------------------------------------------
# 4) RADIATIVE TRANSFER CONFIGURATION
# --------------------------------------------------------------------------------------

Inputs.direct = 0  # 1 = direct-beam, 0 = Diffuse flux
Inputs.aprx_typ = 1  # 1 = Eddington, 2 = Quadrature, 3 = Hemispheric Mean
Inputs.delta = 1  # 1 = Apply delta approximation, 0 = No delta
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
Inputs.add_double = True  # toggle adding-doubling solver

# garder film pour LAPs but too small to be detected
# get varying fresnel layer depth with also undetected depth
# --> how do I vary this depth??...
# then have layer of ice with water until vertical boundary
Inputs.dz = [1e-10, 2e-4, 1]  # thickness of each vertical layer (unit = m) # 0.04
Inputs.nbr_lyr = len(Inputs.dz)  # number of snow layers
Inputs.layer_type = [2, 3, 1]  # 0 = ice grain layer, 1 = solid bubbly ice w/ fresnel layer,  
                               # 2 = liquid water, 3 = solid bubbly ice w/out fresnel layer, 
                               # 4 = solid bubbly ice w/out fresnel layer & water 'bubbles'
Inputs.lwc = 0.06
Inputs.cdom_layer = [0]*len(Inputs.dz)  # Only for layer type == 1
Inputs.rho_layers = [916.999, 780, 780]  # density of each layer (unit = kg m-3)
Inputs.nbr_wvl = 480

# reflectance of underlying surface - set across all wavelengths
Inputs.R_sfc = np.genfromtxt(
    Inputs.dir_base + "Data/OP_data/480band/r_sfc/blue_ice_spectrum_s10290721.csv",
    delimiter="csv",
)
# Inputs.R_sfc = np.array([0 for i in range(Inputs.nbr_wvl)])


# --------------------------------------------------------------------------------------
# 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
# For hexagonal plates or columns of any size choose GeometricOptics
# For sphere, spheroids, koch snowflake with optional water coating choose Mie
# --------------------------------------------------------------------------------------

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
Inputs.grain_rds = [20000,900, 900]  # effective grain or bubble radius
Inputs.rwater = [0]*len(Inputs.dz)  # radius of optional liquid water coating

# For 4:
Inputs.side_length = [1000]*len(Inputs.dz)
Inputs.depth = [1000]*len(Inputs.dz)

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
Inputs.c_factor_GA = 0
Inputs.c_factor_SA = 20

# Set names of files containing the optical properties of these LAPs:
# uncoated BC (Bohren and Huffman, 1983)
Inputs.file_soot1 = "bc_ChCB_rn40_dns1270.nc"
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
Inputs.file_Cook_Greenland_dust_H = "dust_greenland_Cook_CENTRAL_20190911.nc"
# Snow Algae (Cook et al. 2017a, spherical, C nivalis)
Inputs.file_snw_alg = "SA_Chevrollier2022_r8.99.nc"
# Glacier Algae (Cook et al. 2020)
Inputs.file_glacier_algae = "GA_Chevrollier2022_r4.9_L18.8.nc"

# Indicate mass mixing ratios scenarios for each impurity
# default units are ng(species)/g(ice), or ppb
# However, snow and glacier algae can be provided in in cells/mL.
# To use cells/mL for algae, set GA_units == 1.
# The script will loop over the different mixing scenarios
# 1 g/L = 109 ng/L = 106 ng /g

Inputs.mss_cnc_soot1 = [0]* len(Inputs.dz)
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


# --------------------------------------------------------------------------------------
# Error catching: invalid combinations of input variables
# --------------------------------------------------------------------------------------

if Inputs.GA_units == 1 or Inputs.SA_units == 1:  # pylint: disable=W0143
    print("\n****** WARNING ******")
    print(
        "you are expressing at least one of your algal\
    concentrations in cells/mL"
    )
    print(
        " *** This requires your MAC file to be in units of m2/cell***\n\
        our default files are expressed in ppb!"
    )

if Inputs.toon and Inputs.add_double:  # pylint: disable=W0143
    raise ValueError(
        "ERROR: BOTH SOLVERS SELECTED:\
        PLEASE CHOOSE EITHER toon OR add_double"
    )

if Inputs.toon and Inputs.solzen < 40:  # pylint: disable=W0143
    raise ValueError("INVALID SOLZEN: outside valid range for toon solver")

if np.sum(Inputs.layer_type) < 1 and Inputs.add_double:  # pylint: disable=W0143
    # just warn user but let program continue - in some cases
    # AD method preferable (stable over complete range of SZA)
    print("\n****** WARNING ******\n")
    print("No solid ice layers - toon solver is faster")
    print("Toggle toon=True and add_double=False to use it.\n")

if np.sum(Inputs.layer_type) > 0 and Inputs.toon:
    raise ValueError("There are ice layers - use the adding-doubling solver")

if np.sum(Inputs.mss_cnc_snw_alg) != 0:  # pylint: disable=W0143
    # remind user that snow algae optical properties h
    # have not yet been empirically validated
    print("\n****** WARNING ******")
    print("the optical properties for the snow algae are theoretical.")
    print("They were constructed from literature values")
    print("They have not yet been validated empirically.\n\n")

if Inputs.solzen > 89:  # pylint: disable=W0143
    Inputs.solzen = 89
    print(
        "Surface irradiance profiles exist for a solar\
        zenith angle < 90 degrees. Solzen set to 89."
    )
for lyr in range(len(Inputs.dz)):
    if Inputs.cdom_layer[lyr] == 1 and Inputs.layer_type[lyr] == 0: 
       raise ValueError("CDOM option only available for solid ice layers")
# --------------------------------------------------------------------------------------
# CONFIG  AND CALL SNICAR
# --------------------------------------------------------------------------------------

if WRITE_CONFIG_TO_TEXTFILE:
    omitted_fields = [
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
    ]

    with open("./model_config.txt", "w", encoding="utf8") as f:
        for i in Inputs._fields:
            if i not in omitted_fields:
                f.write(i)
                f.write(":  ")
                f.write(str(getattr(Inputs, str(i))))
                f.write("\n")


Outputs = snicar_feeder(Inputs)


# --------------------------------------------------------------------------------------
# PLOTTING AND PRINTING
# --------------------------------------------------------------------------------------
albedo = Outputs.albedo
BBA = Outputs.BBA
wvl = Outputs.wvl
SSA = (
    6
    * (np.array([917] * len(Inputs.dz)) - np.array(Inputs.rho_layers))
    / 917
    / (np.array(Inputs.rho_layers) * np.array(Inputs.grain_rds) * 2 * 10 ** (-6))
)

if SMOOTH:
    from scipy.signal import savgol_filter

    yhat = savgol_filter(albedo, WINDOW_SIZE, POLY_ORDER)
    albedo = yhat

if PRINT_BAND_RATIOS:
    I2DBA = albedo[51] / albedo[46]
    I3DBA = (albedo[46] - albedo[50]) / albedo[55]
    NDCI = ((albedo[50] - albedo[48]) - (albedo[55] - albedo[48])) * (
        (albedo[50] - albedo[48]) / (albedo[55] - albedo[48])
    )
    MCI = (albedo[50] - albedo[46]) / (albedo[50] + albedo[46])
    II = np.log(albedo[36]) / np.log(albedo[66])

    print("\nINDEX VALUES")
    print("2DBA Index: ", I2DBA)
    print("3DBA index: ", I3DBA)
    print("NDCI index: ", NDCI)
    print("MCI index: ", MCI)
    print("Impurity Index: ", II)

if PRINT_BBA:
    print("\nBROADBAND ALBEDO = ", BBA)


rc = {
    "figure.figsize": (8, 6),
    "axes.facecolor": "white",
    "axes.grid": False,
    "grid.color": ".8",
    "xtick.major.width": 1,
    "xtick.major.size": 5,
    "ytick.major.width": 1,
    "ytick.major.size": 5,
    "axes.linewidth": 1.25,
    "font.size": 20,
    "xtick.bottom": True,
    "ytick.left": True,
}
plt.rcParams.update(rc)
plt.plot(wvl[15:230], albedo_test[15:230]) 

plt.plot(wvl[15:230], albedo[15:230], '--') # 0.25m

plt.ylabel("Albedo", fontsize=18), plt.xlabel(
    "Wavelength (µm)", fontsize=18
), plt.xlim(0.3, 2), plt.ylim(0, 1), plt.xticks(fontsize=15), plt.yticks(
    fontsize=15
), 
plt.axvline(
    x=0.95, color="g", linestyle="dashed"
)


if SAVE_FIGS:
    plt.savefig(str(savepath + "spectral_albedo.png"))

