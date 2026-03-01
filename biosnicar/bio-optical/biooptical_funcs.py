#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Joseph Cook, Lou Chevrollier

Contains functions relating to bio-optical model components of BioSNICAR_GO_py.

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from miepython import mie
from scipy.signal import savgol_filter

from biosnicar.classes.bio_optical_config import BioOpticalConfig
from biosnicar.optical_properties.van_diedenhoven import calc_ssa_and_g


def run_biooptical_model(input_file):
    """Executes functions in bio-optical model.

    Calling `run_biooptical_model()` with no arguments runs the full bio-optical model
    with config defined in inputs.yaml. This includes generating new optical properties
    and saving the resuting files to /Data/480band/laps.

    To then incorporate the newly generated impurities into the radiative transfer model
    the new filenames must be used to generate instances of Impurity in setup_snicar().
    This is achieved by adding them to the impurities section of inputs.yaml.

    Args:
        None

    Returns:
        None but saves NetCDFs to /Data/480band/laps.

    """

    bio_optical_config = BioOpticalConfig(input_file)
    abs_cff = get_absorption_cross_section(bio_optical_config)
    k = calculate_k(bio_optical_config, abs_cff)
    (
        wvl_rescaled_BioSNICAR,
        abs_cff_rescaled_BioSNICAR,
        k_rescaled_BioSNICAR,
        n_rescaled_BioSNICAR,
    ) = rescale_480band(bio_optical_config, abs_cff, k)
    plot_k_n_abs_cff(bio_optical_config, abs_cff, k)

    # # --------------------------------------------------------------------------------------
    # # CALCULATIONS OF SCATTERING PROPERTIES
    # # --------------------------------------------------------------------------------------

    assym, ss_alb = calculate_ssps(
        bio_optical_config,
        k_rescaled_BioSNICAR,
        wvl_rescaled_BioSNICAR,
        n_rescaled_BioSNICAR,
    )

    # # --------------------------------------------------------------------------------------
    # # SAVING DATA IN NETCDF
    # # --------------------------------------------------------------------------------------

    net_cdf_updater(
        bio_optical_config,
        assym,
        ss_alb,
        abs_cff_rescaled_BioSNICAR,
        wvl_rescaled_BioSNICAR,
    )

    return


def get_absorption_cross_section(bio_optical_config):

    """Calculates or loads cellular absorption cross section.

    The user can choose whether to calculate the absorption cross section
    for the algal cells by applying the mixing model from Cook et al 2017
    to absorption data for individual pigments (from Dauchet et al 2015),
    or alternatively load a pre-measured/pre-calculated absorption spectrum
    from an external file. The mixing model method has been used in previous
    literature such as Cook et al 2017, 2020, but the default for BioSNICAR
    from v2.0 onwards is to load the absorption spectra from an external file
    provided in this repository, since we were able to make direct empirical
    measurements of the absorption coefficient of intact algal cells.

    Args:
        bio_optical_config: instance of BioOpticalConfig class

    Returns:
        abs_cff: absorption cross section in m2/cell, m2/um3 or m2/mg

    """

    #################
    ## Initialization
    #################

    abs_cff_pigments = pd.DataFrame(
        index=bio_optical_config.wvl * 1000
    )  # storing pigment MACs

    if bio_optical_config.abs_cff_calculated:
        print("abs_cff reconstructed from pigments")
        # open mass absorption coefficients (m2/mg) for each algal pigment
        # from a dictionary.key is pigment name, value is abs coeff in m2/mg
        abs_coeff = 0
        for key, value in bio_optical_config.pigment_data.items():
            abs_pigm = np.array(pd.read_csv(key, header=None)).flatten()  # m2/mg
            abs_cff_pigments[
                str(key.split(bio_optical_config.pigment_dir, 1)[1])[0:-4]
            ] = abs_pigm
            conc = value  # intracellular conc in ng/µm3, ng/cell, or ng/mg
            abs_coeff = abs_coeff + conc * abs_pigm / 1000000  # m2/µm3,m2/cell,m2/mg
        abs_cff = abs_coeff

    elif bio_optical_config.abs_cff_loaded_reconstructed:
        print("abs_cff reconstructed directly loaded")
        abs_cff = np.loadtxt(bio_optical_config.abs_cff_file)  # m2/mg, um3 or cell
        if bio_optical_config.packaging_correction_SA:  # ! applies only from 300nm
            pckg_SA = np.loadtxt(bio_optical_config.dir_pckg + "pckg_SA.csv")
            abs_cff = abs_cff * pckg_SA
        if bio_optical_config.packaging_correction_GA:  # ! applies from 300nm
            pckg_GA = np.loadtxt(bio_optical_config.dir_pckg + "pckg_GA.csv")
            abs_cff = abs_cff * pckg_GA

    elif bio_optical_config.abs_cff_loaded_invivo:
        print("abs_cff in vivo directly loaded")
        abs_cff = np.loadtxt(bio_optical_config.abs_cff_file)  # m2/mg, um3 or cell

    return abs_cff


def calculate_k(bio_optical_config, abs_cff):

    """Calculates imaginary part of refractive index for algal cell.

    Uses the absorption coefficient from file or calculated using mixing
    model to calculate the imaginary part of the refractive index (k) for the
    algal cells. Outside of the visible wavelengths the cells are asumed to
    have k equal to that of water.


    Args:
        bio_optical_config: instance of BioOpticalConfig class
        abs_cff: absorption cross section in m2/cell, m2/um3 or m2/mg

    Returns:
        k: cellular imaginary part of refractive index (unitless)

    """

    k_water_alg = np.loadtxt(bio_optical_config.k_water_dir)
    k_water_alg[0:600] = 0  # accounted for in abs_cff
    xw = (
        0.59 * bio_optical_config.wet_density / 1000
    )  # conversion mass to volume water fraction

    if bio_optical_config.unit == 0:  # abs_cff in m2 / cell
        # units: abs_cff (m2/cell to µm2/cell) / cell volume (um3/cell) * wvl (µm)
        k = (
            xw * k_water_alg
            + abs_cff
            * 10 ** (12)
            * bio_optical_config.wvl
            / (np.pi * 4)
            / bio_optical_config.cell_vol
        )

    if bio_optical_config.unit == 1:  # abs_cff in m2 / µm3
        # units: abs_cff (m2/µm3 to µm2/µm3) * wvl (µm)
        k = xw * k_water_alg + abs_cff * 10 ** (12) * bio_optical_config.wvl / (
            np.pi * 4
        )

    if bio_optical_config.unit == 2:  # abs_cff in m2 / dry mg
        # units: abs_cff (m2/mg to m2/kg) * density (kg m3) * wvl (µm to m)
        k = (
            xw * k_water_alg
            + abs_cff
            * bio_optical_config.dry_density
            * bio_optical_config.wvl
            / (np.pi * 4)
        )

    return k


def rescale_480band(bio_optical_config, abs_cff, k):
    """Rescales abs_cff and refractive index to BioSNICAR resolution.

    The files laoded into the biooptical model have resolution 1nm. This
    function ensures the abs_cff, k and n all have resolution and wavelength
    range equal to that of the BioSNICAR radiative transfer scheme (i.e
    0.205 - 4.995 um in steps on 10 nm)

    Args:
        bio_optical_config: instance of BioOpticalConfig class
        abs_cff: absorption cross section in m2/cell, m2/um3 or m2/mg
        k: cellular imaginary part of refractive index (unitless)

    Returns:
        wvl_rescaled: wavelengths 480 bands
        abs_cff_rescaled: cellular absorption cross section 480 bands (unitless)
        k_rescaled: cellular imaginary part of refractive index 480 bands (unitless)
        n_rescaled: cellular real part of refractive index 480 bands (unitless)

    """

    # rescaling variables to BioSNICAR resolution (10nm)
    wvl_rescaled = bio_optical_config.wvl[5::10]
    abs_cff_rescaled = abs_cff[5::10]
    k_rescaled = k[5::10]
    n_rescaled = (
        bio_optical_config.n_algae * np.ones(np.size(bio_optical_config.wvl))
    )[5::10]

    return wvl_rescaled, abs_cff_rescaled, k_rescaled, n_rescaled


def calculate_ssps(bio_optical_config, k_rescaled, wvl_rescaled, n_rescaled):
    """Calculates single scattering optical properties using Mie or GO scattering codes.

    The user can toggle between using two methods to calculate the single scattering optical
    properties from the real and imaginary parts of the refractive index along with the
    cell size distribution. The BioSNICAR default is to assume monodispersions, but a particle
    size distribution can easily be calculated by running thie bio-optical model over a range
    of particle sizes and calculating an appropriate average. The two schemes availabe are Mie
    scatteriung or geometrical optics. The former assumes spehericalcells and is slow to converge
    for large cells. The latter is fast but assumes a circular cylindrical shape.

    The main calculations used for the geometrical optics mode are based upon the equations of
    Diedenhoven et al (2014), who provided a python script as supplementary material for their
    paper, available at:
    https://www.researchgate.net/publication/259821840_ice_OP_parameterization

    In Mie mode, the optical properties are calculated using Mie scattering
    using Scott Prahl's miepython package https://github.com/scottprahl/miepython.
    This function also plots and saves single scattering properties if toggled.

    Args:
        bio_optical_config: instance of BioOpticalConfig class
        wvl_rescaled: wavelengths 480 bands
        abs_cff_rescaled: cellular absorption cross section 480 bands (unitless)
        k_rescaled: cellular imaginary part of refractive index 480 bands (unitless)
        n_rescaled: cellular real part of refractive index 480 bands (unitless)

    Returns:
        assym: assymetry parameter

        ss_alb: single scattering albedo
    """

    r = bio_optical_config.radius
    L = bio_optical_config.length
    wvl = wvl_rescaled
    n_algae = n_rescaled
    k_algae = k_rescaled

    if bio_optical_config.GO:

        # calculate cylinder dimensions
        diameter = 2 * r
        V = L * (np.pi * r**2)  # volume of cylinder in µm3
        Reff = (V / ((4 / 3) * np.pi)) ** 1 / 3  # effective radius
        Area_total = 2 * (np.pi * r**2) + (2 * np.pi * r) * L
        Area = Area_total / 4  # projected area
        ar = diameter / L  # aspect ratio

        # SSA and asymmetry parameter via shared Van Diedenhoven parameterization
        ss_alb, assym = calc_ssa_and_g(ar, V, Area, n_algae, k_algae, wvl)

        # Absorption cross section and size parameter (for optional plots)
        absXS = (1 - np.exp(-4 * np.pi * k_algae * V / Area * wvl)) * Area
        X = (2 * np.pi * Reff) / wvl

        if bio_optical_config.plot_ssps:

            plt.figure(1)
            plt.plot(wvl, ss_alb, label="{}x{}".format(r, L)),
            plt.xlim(0.3, 1.4), plt.ylabel("SSA"),
            plt.xlabel("Wavelength (um)"), plt.grid(False),
            plt.legend(loc="best", ncol=2)

            plt.figure(2)
            plt.plot(wvl, assym, label="{}x{}".format(r, L)),
            plt.xlim(0.3, 1.4), plt.ylabel("Assymetry Parameter"),
            plt.xlabel("Wavelength (um)"), plt.grid(False),
            plt.legend(loc="best", ncol=2)

            plt.figure(3)
            plt.plot(wvl, absXS, label="{}x{}".format(r, L)),
            plt.xlim(0.3, 1.4), plt.ylabel("Absorption Cross Section"),
            plt.xlabel("Wavelength (um)"), plt.grid(False),
            plt.legend(loc="best", ncol=2)

            plt.figure(4)
            plt.plot(wvl, X, label="{}x{}".format(r, L)),
            plt.ylabel("Size Parameter X"), plt.xlabel("Wavelength (um)"),
            plt.grid(False), plt.legend(loc="best", ncol=2)

            plt.show()

            if bio_optical_config.savefig_ssps:
                plt.savefig(
                    str(bio_optical_config.savepath + "SSA_{}x{}.jpg".format(r, L))
                )
                plt.savefig(
                    str(
                        bio_optical_config.savepath
                        + "AssymetryParam_{}x{}.jpg".format(r, L)
                    )
                )
                plt.savefig(
                    str(bio_optical_config.savepath + "AbsXS_{}x{}.jpg".format(r, L))
                )
                plt.savefig(
                    str(
                        bio_optical_config.savepath
                        + "X_SizeParam_{}x{}.jpg".format(r, L)
                    )
                )

        if bio_optical_config.report_dims:
            print("cell diameter = ", np.round(diameter, 2), " (micron)")
            print("cell length = ", L, " (micron)")
            print("aspect ratio = ", ar)
            print("cell volume = ", np.round(V, 2), " (micron^3)")
            print("Effective radius = ", np.round(Reff, 2), " (micron)")
            print("projected area = ", np.round(Area, 2))
            print()  # line break

    if bio_optical_config.Mie:

        X = 2 * np.pi * r / wvl  # unitless
        qext, qsca, qback, g = mie(n_algae - 1j * k_algae, X)
        qabs = qext - qsca
        assym = g
        ss_alb = qsca / qext

        if bio_optical_config.plot_ssps:
            qqabs = qabs * np.pi * r**2  # calculate cross section from efficiency
            qqsca = qsca * np.pi * r**2  # calculate cross section from efficiency
            qqext = qext * np.pi * r**2  # calculate cross section from efficiency
            plt.figure(1)
            plt.plot(wvl, qqabs, "b", label="absorption cross section")
            plt.plot(wvl, qqsca, "g", label="scattering cross section")
            plt.plot(wvl, qqext, "r", label="extinction cross section")
            plt.xlim(0.2, 2.5)

            plt.ylabel(r"Cross Section ($\mu$m$^2$)")
            plt.xlabel("Wavelength (nm)")
            plt.legend(loc="best")
            plt.show()

            if bio_optical_config.savefig_ssps:
                plt.savefig(
                    str(bio_optical_config.savepath + "Mie_params.jpg"), dpi=150
                )
                plt.show()

    return assym, ss_alb


def plot_k_n_abs_cff(bio_optical_config, abs_cff, k):

    """Optionally plots and saves figures and files for abs_cff, n and k.

    Args:
        bio_optical_config: instance of BioOpticalConfig class
        abs_cff: absorption cross section in m2/cell, m2/um3 or m2/mg
        k: cellular imaginary part of refractive index (unitless)

    Returns:
        None
    """

    if bio_optical_config.smooth:
        yhat = savgol_filter(
            abs_cff, bio_optical_config.window_size, bio_optical_config.poly_order
        )  # window size 51, polynomial order 3
        abs_cff = yhat

    # optionally save files to savepath
    if bio_optical_config.savefiles_n_k_abs_cff:  # optional save dataframe to csv files
        pd.DataFrame(k).to_csv(
            str(bio_optical_config.savepath + "k.csv"), header=None, index=False
        )
        pd.DataFrame(abs_cff).to_csv(
            str(bio_optical_config.savepath + "abs_cff.csv"), header=None, index=False
        )
        pd.DataFrame(
            (bio_optical_config.n_algae * np.ones(np.size(bio_optical_config.wvl)))
        ).to_csv(str(bio_optical_config.savepath + "n.csv"), header=None, index=False)

    # optionally plot figures to interative window
    if bio_optical_config.plot_k_abs_cff:
        plt.figure(figsize=(7, 9))
        plt.subplot(2, 1, 1)
        plt.plot(bio_optical_config.wvl[100:600] * 1000, abs_cff[100:600])
        plt.xticks(fontsize=15), plt.yticks(fontsize=15)
        plt.xlabel("Wavelength (nm)", fontsize=15),
        plt.ylabel(
            str("abs_cff (m$2$ cell$^{-1}$, m$^2$ µm$3$ or m$^2$ ng$3$ ))"), fontsize=15
        )
        plt.tight_layout()

        plt.subplot(2, 1, 2)
        plt.plot(bio_optical_config.wvl[100:600] * 1000, k[100:600])
        plt.xticks(fontsize=15), plt.yticks(fontsize=15)
        plt.xlabel("Wavelength (nm)", fontsize=15), plt.ylabel("k", fontsize=15)
        plt.tight_layout()
        plt.show()

        # optionally save plots to savepath
        if bio_optical_config.saveplots_k_abs_cff:
            plt.show()
            plt.savefig(str(bio_optical_config.savepath + "/abs_cffandK.png"))


def net_cdf_updater(bio_optical_config, assym, ss_alb, abs_cff_rescaled, wvl_rescaled):
    """Optionally saves netcdf file with optical properties usable in BioSNICAR.

    Saves the calculate optical properties to a NetCDF file with formatting that allows the
    data to be loaded into BioSNICAR.

    Args:
        bio_optical_config: instance of BioOpticalConfig class
        abs_cff_rescaled: absorption cross section in m2/cell, m2/um3 or m2/mg
        assym: assym parameter
        ss_alb: single scattering albedo
        wvl_rescaled: wavelengths

    Returns:
        None
    """

    if bio_optical_config.save_netcdf:
        algfile = pd.DataFrame()
        algfile["asm_prm"] = np.squeeze(assym)
        algfile["ss_alb"] = np.squeeze(ss_alb)
        algfile["ext_cff_mss"] = abs_cff_rescaled
        algfile = algfile.to_xarray()
        algfile.attrs["medium_type"] = "air"
        if bio_optical_config.GO:
            algfile.attrs["description"] = (
                "Optical properties for glacier algal cell: cylinder of radius "
                "{}um and length {}um".format(
                    str(bio_optical_config.radius), str(bio_optical_config.length)
                )
            )
            algfile.attrs["side_length_um"] = bio_optical_config.length
        if bio_optical_config.Mie:
            algfile.attrs["description"] = (
                "Optical properties for snow algal "
                "cell: sphere of radius {}um".format(str(bio_optical_config.radius))
            )
        algfile.attrs["psd"] = "monodisperse"
        algfile.attrs["density_kg_m3"] = bio_optical_config.wet_density
        algfile.attrs["wvl"] = wvl_rescaled
        algfile.attrs["information"] = bio_optical_config.information
        algfile.to_netcdf(
            str(
                bio_optical_config.savepath_netcdf
                + bio_optical_config.filename_netcdf
                + ".nc"
            ),
            mode="w",
        )

    return


if __name__ == "__main__":
    pass
