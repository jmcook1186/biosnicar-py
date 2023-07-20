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

    ssa_list = []
    g_list = []
    abs_xs_list = []
    mac_list = []

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
    delta = 0.3

    for i in np.arange(0, len(wavelengths), 1):
        mr = reals[i]
        mi = imags[i]
        wl = wavelengths[i]

        # ------------------------------------------------
        # ---------- input tables (see Figs. 4 and 7) ----
        # ------------------------------------------------
        # SSA parameterization
        a = [0.457593, 20.9738]  # for ar=1

        # SSA correction for AR != 1 (Table 2)
        nc1 = 3
        nc2 = 4
        c_ij = np.zeros(nc1 * nc2 * 2).reshape((nc1, nc2, 2))
        #   ---------- Plates ----------
        c_ij[:, 0, 0] = [0.000527060, 0.309748, -2.58028]
        c_ij[:, 1, 0] = [0.00867596, -0.650188, -1.34949]
        c_ij[:, 2, 0] = [0.0382627, -0.198214, -0.674495]
        c_ij[:, 3, 0] = [0.0108558, -0.0356019, -0.141318]
        #   --------- Columns ----------
        c_ij[:, 0, 1] = [0.000125752, 0.387729, -2.38400]
        c_ij[:, 1, 1] = [0.00797282, 0.456133, 1.29446]
        c_ij[:, 2, 1] = [0.00122800, -0.137621, -1.05868]
        c_ij[:, 3, 1] = [0.000212673, 0.0364655, 0.339646]

        # diffraction g parameterization
        b_gdiffr = [-0.822315, -1.20125, 0.996653]

        # raytracing g parameterization ar=1
        p_a_eq_1 = [0.780550, 0.00510997, -0.0878268, 0.111549, -0.282453]

        # ---- g correction for AR != 1 (Also applied to AR=1 as plate) (Table 3)
        nq1 = 3
        nq2 = 7
        q_ij = np.zeros(nq1 * nq2 * 2).reshape((nq1, nq2, 2))
        #   ---------- Plates ----------
        q_ij[:, 0, 0] = [-0.00133106, -0.000782076, 0.00205422]
        q_ij[:, 1, 0] = [0.0408343, -0.00162734, 0.0240927]
        q_ij[:, 2, 0] = [0.525289, 0.418336, -0.818352]
        q_ij[:, 3, 0] = [0.443151, 1.53726, -2.40399]
        q_ij[:, 4, 0] = [0.00852515, 1.88625, -2.64651]
        q_ij[:, 5, 0] = [-0.123100, 0.983854, -1.29188]
        q_ij[:, 6, 0] = [-0.0376917, 0.187708, -0.235359]
        #   ---------- Columns ----------
        q_ij[:, 0, 1] = [-0.00189096, 0.000637430, 0.00157383]
        q_ij[:, 1, 1] = [0.00981029, 0.0409220, 0.00908004]
        q_ij[:, 2, 1] = [0.732647, 0.0539796, -0.665773]
        q_ij[:, 3, 1] = [-1.59927, -0.500870, 1.86375]
        q_ij[:, 4, 1] = [1.54047, 0.692547, -2.05390]
        q_ij[:, 5, 1] = [-0.707187, -0.374173, 1.01287]
        q_ij[:, 6, 1] = [0.125276, 0.0721572, -0.186466]

        # --------- refractive index correction of asymmetry parameter
        c_g = np.zeros(4).reshape(2, 2)
        c_g[:, 0] = [0.96025050, 0.42918060]
        c_g[:, 1] = [0.94179149, -0.21600979]
        # ---- correction for absorption
        s = [1.00014, 0.666094, -0.535922, -11.7454, 72.3600, -109.940]
        u = [-0.213038, 0.204016]

        # -------- selector for plates or columns
        if ar > 1.0:
            col_pla = 1  # columns
        else:
            col_pla = 0  # plates & compacts

        # ------------------------------------------------
        # ------------ Size parameters -------------------
        # ------------------------------------------------

        # --- absorption size parameter (Fig. 4, box 1)
        Chi_abs = mi / wl * V / Area

        # ----- scattering size parameter (Fig. 7, box 1)
        Chi_scat = 2.0 * np.pi * np.sqrt(Area / np.pi) / wl

        # ------------------------------------------------
        # ------------ SINGLE SCATTERING ALBEDO ----------
        # ------------------------------------------------

        if Chi_abs > 0:
            w_1 = 1.0 - a[0] * (
                1.0 - np.exp(-Chi_abs * a[1])
            )  # for AR=1 (Fig. 4, box 2)
            l = np.zeros(nc1)
            for i in range(nc2):
                l[:] += c_ij[:, i, col_pla] * np.log10(ar) ** i  # (Fig. 4, box 3)
            D_w = (
                l[0]
                * np.exp(-((np.log(Chi_abs) - l[2]) ** 2) / (2.0 * l[1] ** 2))
                / (Chi_abs * l[1] * np.sqrt(2.0 * np.pi))
            )  # (Fig. 4, box 3)
            w = w_1 + D_w  # (Fig. 4, box 4)
        else:
            w = 1.0

        # ------------------------------------------------
        # --------------- ASYMMETRY PARAMETER ------------
        # ------------------------------------------------

        # diffraction g
        g_diffr = (
            b_gdiffr[0] * np.exp(b_gdiffr[1] * np.log(Chi_scat)) + b_gdiffr[2]
        )  # (Fig. 7, box 2)
        g_diffr = max([g_diffr, 0.5])

        # raytracing g at 862 nm
        g_1 = 0.0
        for i in range(len(p_a_eq_1)):
            g_1 += p_a_eq_1[i] * delta**i  # (Fig. 7, box 3)

        p_delta = np.zeros(nq1)
        for i in range(nq2):
            p_delta += q_ij[:, i, col_pla] * np.log10(ar) ** i  # (Fig. 7, box 4)
        Dg = 0.0
        for i in range(nq1):
            Dg += p_delta[i] * delta**i  # (Fig. 7, box 4)
        g_rt = 2.0 * (g_1 + Dg) - 1.0  # (Fig. 7, box 5)

        # --------- refractive index correction of asymmetry parameter (Fig. 7, box 6)
        epsilon = c_g[0, col_pla] + c_g[1, col_pla] * np.log10(ar)
        mr1 = 1.3038  # reference value @ 862 nm band
        C_m = abs(
            (mr1 - epsilon) / (mr1 + epsilon) * (mr + epsilon) / (mr - epsilon)
        )  # abs function added according to corrigendum to the original paper

        # ---- correction for absorption (Fig. 7, box 7)
        if Chi_abs > 0:
            C_w0 = 0.0
            for i in range(len(s)):
                C_w0 += s[i] * (1.0 - w) ** i
            k = np.log10(ar) * u[col_pla]
            C_w1 = k * w - k + 1.0
            C_w = C_w0 * C_w1
        else:
            C_w = 1.0

        # raytracing g at required wavelength
        g_rt_corr = g_rt * C_m * C_w  # (Fig. 7, box 9)

        # ----- Calculate total asymmetry parameter and check g_tot <= 1 (Fig. 7, box 9)
        g_tot = 1.0 / (2.0 * w) * ((2.0 * w - 1.0) * g_rt_corr + g_diffr)
        g_tot = min([g_tot, 1.0])

        absXS = Area * (1 - ((np.exp(-4 * np.pi * mi * V)) / (Area * wl)))
        MAC = (
            absXS / V * 914
        )  # divide by volume*mass to give mass absorption coefficient

        ssa_list.append(w)
        g_list.append(g_tot)
        abs_xs_list.append(absXS)
        mac_list.append(MAC)

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
