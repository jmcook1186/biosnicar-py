#!/usr/bin/python

"""Calculates optical properties of large hexagonal ice grains.

This file calculates the optial properties (single scattering albedo, assymetry
parameter, mass absorption coefficient and extinction, scattering and absorption cross
sections) for ice grains shaped as arbitrarily large hexagonal plates or columns.
The optical propertiesare then saved into netCDF files in the correct format for
loading into BioSNICAR.

The main function calc_optical_params() is based upon the equations of
Diedenhoven et al (2014) who provided a python script as supplementary material
for their paper. The original code can be downloaded from:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

The optical properties are calculated using a parameterization of geometric optics
calculations (Macke et al., JAS, 1996).


There are no user defined inputs for the preprocessing function, it can simply be
run as

reals, imags, wavelengths = preprocess()

The calc_optical_params() fnction takes several inputs. reals, imags and wavelengths
are output by preprocess() and side_length and depth are user defined. These are the two
parameters that control the dimensions of the ice crystals. Side_length is the length
in microns of one side of the hexagnal face of the crystal, depth is the column length
also in microns describing the z dimension. The code then calculates volume, apotherm,
aspect ratio, area etc inside the function. The optical parameters are returned.
Optional plots and printed values for the optical params are provided by setting
plots to true and the dimensions of the crystals can be reported by setting
report_dims to true in the function call.

The final function, net_cdf_updater() is used to dump the optical parameters and
metadata into a netcdf file and save it into the working directory to be used as
a lookup library for the two-stream radiative transfer model BoSNICAR_GO.

The function calls are provided at the bottom of this script in a loop, where the
user can define the range of side lengths and depths to be looped over.

NOTE: The extinction coefficient in the current implementation is 2 for all size
parameters as assumed in the conventional geometric optics approximation.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from biosnicar.optical_properties.van_diedenhoven import calc_ssa_and_g

# Set paths
SAVEPATH = "./Data/GO_files/480band/"
DATAPATH = "./Data/rfidx_ice.nc"
RI_SOURCE = 2


def preprocess_RI(ri_source, path_to_ri):
    """Preprocessing of wavelength and RI data.

    Preprocessing function that ensures the wavelengths and real/imaginary
    parts of the refractive index for ice is provided in the correct waveband and correct
    spectral resolution to interface with BioSNICAR. The refractive indices are taken
    from Warren and Brandt 2008.

    Grabs appropriates wavelengths, real and imaginary parts of ice
    refractive index. The source of the refractive index data is
    controlled by var "ri_source" where 0 = Warren 1984, 1 = Warren 2008
    and 2 = Picard 2016.

    These are then passed as numpy arrays to the Geometrical Optics function.

    Args:
        ri_source: choice of refractive index
        path_to_ri: path to directory containing RI data

    Returns:
        reals: numpy array of real parts of RI by wavelength
        imags: numpy array of imaginary parts of RI by wavelength
        wavelengths: numpy array of wavelengths (um)
    """

    refidx = xr.open_dataset(path_to_ri)
    wavelengths = refidx["wvl"].values

    if ri_source == 0:
        reals = refidx["re_Wrn84"].values
        imags = refidx["im_Wrn84"].values

    elif ri_source == 1:
        reals = refidx["re_Wrn08"].values
        imags = refidx["im_Wrn08"].values

    elif ri_source == 2:
        reals = refidx["re_Pic16"].values
        imags = refidx["im_Pic16"].values

    return reals, imags, wavelengths


def calc_optical_params(
    side_length,
    depth,
    reals,
    imags,
    wavelengths,
    plots=False,
    report_dims=False,
):
    """Calculates single scattering optical properties.

    Van Diedenhoven's parameterisation is used to calculate
    the single scatterign optical properties of hexagonal
    ice columns of given dimensions.

    Args:
        side_length: length of side of hexagonal face (um)
        depth: length of hexagonal column (um)
        reals: numpy array of real parts of RI by wavelength
        imags: numpy array of imaginary parts of RI by wavelength
        wavelengths: numpy array of wavelenmgths (um)
        plots: Boolean to toggle plotting OPs
        report_dims: Boolean to toggle printing OP data to terminal

    Returns:
        g_list: assymetry parameter
        ssa_list: single scattering albedo
        mac_list: mass absorption coefficient
        depth: length of hexagional column (um)
        side_length: length of side of hexagonal face (um)
        diameter: diameter across hexaginal face.

    """

    V = 1.5 * np.sqrt(3) * side_length**2 * depth  # volume
    Area_total = (
        3 * side_length * (np.sqrt(3) * side_length + depth * 2)
    )  # total surface area
    Area = Area_total / 4  # projected area
    apothem = (2 * Area) / (
        depth * 6
    )  # apothem is distance from centre point to midpoint of a side for hexagon
    diameter = 2 * apothem  # midpoint of one side to midpoint of opposite side

    ar = depth / side_length

    # SSA and asymmetry parameter via shared Van Diedenhoven parameterization
    ssa_list, g_list = calc_ssa_and_g(ar, V, Area, reals, imags, wavelengths)

    # Absorption cross section and mass absorption coefficient
    absXS = Area * (1 - np.exp(-4 * np.pi * imags * V) / (Area * wavelengths))
    mac_list = absXS / V * 914

    if plots:
        plt.figure(1)
        plt.plot(wavelengths, ssa_list), plt.ylabel("SSA"), plt.xlabel(
            "Wavelength (um)"
        ), plt.grid(b=None)
        plt.figure(2)
        plt.plot(wavelengths, g_list), plt.ylabel("Assymetry Parameter"), plt.xlabel(
            "Wavelength (um)"
        ), plt.grid(b=None)
        plt.figure(3)
        plt.plot(wavelengths, mac_list), plt.ylabel(
            "Mass Absorption Cross Section"
        ), plt.xlabel("Wavelength (um)"), plt.grid(b=None)

    if report_dims:
        print("Width of hexagonal plane = ", np.round(diameter / 10000, 2), " (cm)")
        print("depth of hexagonal column = ", depth / 10000, " (cm)")
        print("aspect ratio = ", ar)
        print("ice crystal volume = ", np.round(V * 1e-12, 2), " (cm^3)")

    return g_list, ssa_list, mac_list, depth, side_length, diameter


def net_cdf_updater(
    ri_source, savepath, g_list, ssa_list, mac_list, depth, side_length, density
):
    """Updates a template NetCDF file with new OP data.

    Args:
        ri_source: chocie of refractive index file
        savepath: path to save output data
        g_list: asymmetry parameter
        ssa_list: single scattering albedo
        mac_list: mass absorption coefficient
        depth: length of hexagional column (um)
        side_length: length of side of hexagonal face (um)
        density: density of material in kg/m3.

    Returns:
        None but saves NetCDF file to savepath

    """

    filepath_in = savepath
    mac_in = np.squeeze(mac_list)
    ssa_in = np.squeeze(ssa_list)
    g_in = np.squeeze(g_list)

    if ri_source == 0:
        stb1 = "ice_Wrn84/"
        stb2 = "ice_Wrn84_"

    elif ri_source == 1:
        stb1 = "ice_Wrn08/"
        stb2 = "ice_Wrn08_"

    elif ri_source == 2:
        stb1 = "ice_Pic16/"
        stb2 = "ice_Pic16_"

    icefile = pd.DataFrame()
    icefile["asm_prm"] = g_in
    icefile["ss_alb"] = ssa_in
    icefile["ext_cff_mss"] = mac_in
    icefile = icefile.to_xarray()
    icefile.attrs["medium_type"] = "air"
    icefile.attrs[
        "description"
    ] = f"""Optical properties for ice grain: hexagonal column of side
    length {side_length}um and length {depth}um"""
    icefile.attrs["psd"] = "monodisperse"
    icefile.attrs["side_length_um"] = depth
    icefile.attrs["density_kg_m3"] = density
    icefile.attrs[
        "origin"
    ] = "Optical properties derived from geometrical optics calculations"
    icefile.to_netcdf(
        str(filepath_in + stb1 + stb2 + "{}_{}.nc".format(str(side_length), str(depth)))
    )

    return


# --------------------------------------------------------------------------------------
# FUNCTON CALLS
# --------------------------------------------------------------------------------------

# reals, imags, wavelengths = preprocess_RI(RI_SOURCE, DATAPATH)

# for side_length in np.arange(2000, 11000, 1000):
#     for depth in np.arange(2000, 31000, 1000):

#         (
#             g_list,
#             ssa_list,
#             mac_list,
#             depth,
#             side_length,
#             diameter,
#         ) = calc_optical_params(
#             side_length, depth, reals, imags, wavelengths, plots=False, report_dims=True
#         )

#         net_cdf_updater(
#             RI_SOURCE, SAVEPATH, g_list, ssa_list, mac_list, depth, side_length, 917
#         )


if __name__ == "__main__":
    pass
