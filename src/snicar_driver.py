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
from classes import *

with open("/home/joe/Code/BioSNICAR_GO_PY/src/impurity_config.yaml", "r") as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

impurities = []

for i, id in enumerate(cfg["impurities"]):
    name = cfg["impurities"][id]["name"]
    file = cfg["impurities"][id]["file"]
    cfactor = cfg["impurities"][id]["cfactor"]
    coated = cfg["impurities"][id]["coated"]
    unit = cfg["impurities"][id]["unit"]
    conc = cfg["impurities"][id]["conc"]
    impurities.append(Impurity(dir_base, file, coated, cfactor, unit, name, conc))

ice = Ice(dir_base)
rt = RTConfig()
illumination = illumination()

illumination.calculate_flx_slr()

print(illumination.Fd)








# Inputs = c.namedtuple(
#     "Inputs",
#     [
#         "dir_base",
#         "verbosity",
#         "rf_ice",
#         "incoming_i",
#         "direct",
#         "layer_type",
#         "cdom_layer",
#         "aprx_typ",
#         "delta",
#         "solzen",
#         "toon",
#         "add_double",
#         "R_sfc",
#         "dz",
#         "rho_layers",
#         "grain_rds",
#         "side_length",
#         "depth",
#         "rwater",
#         "nbr_lyr",
#         "nbr_aer",
#         "grain_shp",
#         "shp_fctr",
#         "grain_ar",
#         "GA_units",
#         "SA_units",
#         "c_factor_GA",
#         "c_factor_SA",
#         "mss_cnc_soot1",
#         "mss_cnc_soot2",
#         "mss_cnc_brwnC1",
#         "mss_cnc_brwnC2",
#         "mss_cnc_dust1",
#         "mss_cnc_dust2",
#         "mss_cnc_dust3",
#         "mss_cnc_dust4",
#         "mss_cnc_dust5",
#         "mss_cnc_ash1",
#         "mss_cnc_ash2",
#         "mss_cnc_ash3",
#         "mss_cnc_ash4",
#         "mss_cnc_ash5",
#         "mss_cnc_ash_st_helens",
#         "mss_cnc_Skiles_dust1",
#         "mss_cnc_Skiles_dust2",
#         "mss_cnc_Skiles_dust3",
#         "mss_cnc_Skiles_dust4",
#         "mss_cnc_Skiles_dust5",
#         "mss_cnc_GreenlandCentral1",
#         "mss_cnc_GreenlandCentral2",
#         "mss_cnc_GreenlandCentral3",
#         "mss_cnc_GreenlandCentral4",
#         "mss_cnc_GreenlandCentral5",
#         "mss_cnc_Cook_Greenland_dust_L",
#         "mss_cnc_Cook_Greenland_dust_C",
#         "mss_cnc_Cook_Greenland_dust_H",
#         "mss_cnc_snw_alg",
#         "mss_cnc_glacier_algae",
#         "file_soot1",
#         "file_soot2",
#         "file_brwnC1",
#         "file_brwnC2",
#         "file_dust1",
#         "file_dust2",
#         "file_dust3",
#         "file_dust4",
#         "file_dust5",
#         "file_ash1",
#         "file_ash2",
#         "file_ash3",
#         "file_ash4",
#         "file_ash5",
#         "file_ash_st_helens",
#         "file_Skiles_dust1",
#         "file_Skiles_dust2",
#         "file_Skiles_dust3",
#         "file_Skiles_dust4",
#         "file_Skiles_dust5",
#         "file_GreenlandCentral1",
#         "file_GreenlandCentral2",
#         "file_GreenlandCentral3",
#         "file_GreenlandCentral4",
#         "file_GreenlandCentral5",
#         "file_Cook_Greenland_dust_L",
#         "file_Cook_Greenland_dust_C",
#         "file_Cook_Greenland_dust_H",
#         "file_snw_alg",
#         "file_glacier_algae",
#         "tau",
#         "g",
#         "SSA",
#         "mu_not",
#         "nbr_wvl",
#         "wvl",
#         "Fs",
#         "Fd",
#         "L_snw",
#         "flx_slr",
#     ],
# )



# nbr_aer = len(impurities)

# ## ERROR CHECKING
# # if (
# #     (np.sum(Inputs.mss_cnc_Cook_Greenland_dust_L) > 0)
# #     or (np.sum(Inputs.mss_cnc_Cook_Greenland_dust_C) > 0)
# #     or (np.sum(Inputs.mss_cnc_Cook_Greenland_dust_H) > 0)
# # ):
# #     raise ValueError(
# #         "The Greenland Dusts are not available in this release.\
# #                      They will be included in the next version."
# #     )

# # if Inputs.GA_units == 1 or Inputs.SA_units == 1:  # pylint: disable=W0143
# #     print("\n****** WARNING ******")
# #     print(
# #         "you are expressing at least one of your algal\
# #     concentrations in cells/mL"
# #     )
# #     print(
# #         " *** This requires your MAC file to be in units of m2/cell***\n\
# #         our default files are expressed in ppb!"
# #     )

# # if Inputs.toon and Inputs.add_double:  # pylint: disable=W0143
# #     raise ValueError(
# #         "ERROR: BOTH SOLVERS SELECTED:\
# #         PLEASE CHOOSE EITHER toon OR add_double"
# #     )

# # if Inputs.toon and Inputs.solzen < 40:  # pylint: disable=W0143
# #     raise ValueError("INVALID SOLZEN: outside valid range for toon solver")

# # if np.sum(Inputs.layer_type) < 1 and Inputs.add_double:  # pylint: disable=W0143
# #     # just warn user but let program continue - in some cases
# #     # AD method preferable (stable over complete range of SZA)
# #     print("\n****** WARNING ******\n")
# #     print("No solid ice layers - toon solver is faster")
# #     print("Toggle toon=True and add_double=False to use it.\n")

# # if np.sum(Inputs.layer_type) > 0 and Inputs.toon:
# #     raise ValueError("There are ice layers - use the adding-doubling solver")

# # if np.sum(Inputs.mss_cnc_snw_alg) != 0:  # pylint: disable=W0143
# #     # remind user that snow algae optical properties h
# #     # have not yet been empirically validated
# #     print("\n****** WARNING ******")
# #     print("the optical properties for the snow algae are theoretical.")
# #     print("They were constructed from literature values")
# #     print("They have not yet been validated empirically.\n\n")

# # if Inputs.solzen > 89:  # pylint: disable=W0143
# #     Inputs.solzen = 89
# #     print(
# #         "Surface irradiance profiles exist for a solar\
# #         zenith angle < 90 degrees. Solzen set to 89."
# #     )



# Outputs = snicar_feeder(Inputs)


# # --------------------------------------------------------------------------------------
# # PLOTTING AND PRINTING
# # --------------------------------------------------------------------------------------
# albedo = Outputs.albedo
# BBA = Outputs.BBA
# wvl = Outputs.wvl
# SSA = (
#     6
#     * (np.array([917] * len(Inputs.dz)) - np.array(Inputs.rho_layers))
#     / 917
#     / (np.array(Inputs.rho_layers) * np.array(Inputs.grain_rds) * 2 * 10 ** (-6))
# )

# if SMOOTH:
#     from scipy.signal import savgol_filter

#     yhat = savgol_filter(albedo, WINDOW_SIZE, POLY_ORDER)
#     albedo = yhat

# if PRINT_BAND_RATIOS:
#     I2DBA = albedo[51] / albedo[46]
#     I3DBA = (albedo[46] - albedo[50]) / albedo[55]
#     NDCI = ((albedo[50] - albedo[48]) - (albedo[55] - albedo[48])) * (
#         (albedo[50] - albedo[48]) / (albedo[55] - albedo[48])
#     )
#     MCI = (albedo[50] - albedo[46]) / (albedo[50] + albedo[46])
#     II = np.log(albedo[36]) / np.log(albedo[66])

#     print("\nINDEX VALUES")
#     print("2DBA Index: ", I2DBA)
#     print("3DBA index: ", I3DBA)
#     print("NDCI index: ", NDCI)
#     print("MCI index: ", MCI)
#     print("Impurity Index: ", II)

# if PRINT_BBA:
#     print("\nBROADBAND ALBEDO = ", BBA)

# plt.style.use("seaborn")
# sns.set_style("white")
# rc = {
#     "figure.figsize": (8, 6),
#     "axes.facecolor": "white",
#     "axes.grid": False,
#     "grid.color": ".8",
#     "xtick.major.width": 1,
#     "xtick.major.size": 5,
#     "ytick.major.width": 1,
#     "ytick.major.size": 5,
#     "axes.linewidth": 1.25,
#     "font.size": 20,
#     "xtick.bottom": True,
#     "ytick.left": True,
# }
# plt.rcParams.update(rc)
# plt.plot(wvl, albedo), plt.ylabel("Albedo", fontsize=18), plt.xlabel(
#     "Wavelength (µm)", fontsize=18
# ), plt.xlim(0.3, 2.5), plt.ylim(0, 1), plt.xticks(fontsize=15), plt.yticks(
#     fontsize=15
# ), plt.axvline(
#     x=0.68, color="g", linestyle="dashed"
# )

# if SHOW_FIGS:
#     plt.show()

# if SAVE_FIGS:
#     plt.savefig(str(savepath + "spectral_albedo.png"))
