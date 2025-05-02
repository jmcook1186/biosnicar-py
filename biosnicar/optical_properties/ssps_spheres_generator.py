#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Lou-Anne Chevrollier, University of Aarhus

This code is used to generate the optical properties of an air bubble or grain
size distribution of effective radius re and save them in a netcdf file that
can directly be used in the SNICAR/biosnicar model. The methodology follows that of 
Flanner et al. 2021 (https://doi.org/10.5194/gmd-14-7673-2021), except that 
the Mie solver used is that of Wiscombe 1979 (doi:10.5065/D6ZP4414),
implemented in python by Scott Prahl in the package miepython, available at: 
https://miepython.readthedocs.io/en/latest/.

Temporary files are first created in a temporary folder to store the OPs of 
individual spheres and then mixed to be finally saved in the data folder of 
the SNICAR model.

"""

import numpy as np
import pandas as pd
import miepython as mie
import xarray as xr
import glob
import time

##############################################################################
# Inputs of the Mie solver
##############################################################################

def set_inputs_to_mie_solver():
    wvl = np.arange(0.205e-6, 5e-6, 0.01e-6)
    ref_index_ice = xr.open_dataset("Data/OP_data/480band/rfidx_ice.nc")
    ref_index_water = pd.read_csv(
        "Data/OP_data/480band/refractive_index_water_273K_Rowe2020.csv"
    )
    n_water = ref_index_water.n.values
    k_water = ref_index_water.k.values
    n_ice = ref_index_ice.re_Pic16.values
    k_ice = ref_index_ice.im_Pic16.values
    n_air = np.ones(480) + 1e-6 * (
        (77.46 + 0.459 / ((wvl * 1e6) ** 2)) * 101325 * 0.01 / 273.16
    )
    description = "Wiscombe 1979 solver implemented by Scott Prahl (miepython)"


    path_to_save_temp = "Data/OP_data/480band/tmp8"  # path to temporary individual ops
    air_bbl, water_bbl, ice_grain, water_grain = True, False, False, False
    small_sizes = True  # small is below 5000um for bubbles and 1500um for grains

    if air_bbl:
        path_to_save_ops = "Data/OP_data/480band/bubbly_ice_files"
        filename = "bbl"
        medium_type = "ice_Pic16"
        particle_type = "air_stp"
        particle_density = 1.293
        n_in = n_air
        k_in = np.zeros(480)
        n_ext = n_ice
        if small_sizes:
            sz_min = 5e-08
            sz_max = 3e-02
            sz_nbr = 10000
        else:
            sz_min = 5e-05
            sz_max = 1e-01
            sz_nbr = 500

    if water_bbl:
        path_to_save_ops = "Data/OP_data/480band/bubbly_ice_files"
        filename = "bbl_water"
        medium_type = "ice_Pic16"
        particle_type = "bbl_water"
        particle_density = 1000
        n_in = n_water
        k_in = k_water
        n_ext = n_ice
        if small_sizes:
            sz_min = 5e-08
            sz_max = 3e-02
            sz_nbr = 10000
        else:
            sz_min = 5e-05
            sz_max = 1e-01
            sz_nbr = 500

    if ice_grain:
        path_to_save_ops = "../Data/OP_data/480band/ice_spherical_grains/ice_Pic16"
        filename = "ice_Pic16"
        medium_type = "air_stp"
        particle_type = "ice_Pic16"
        particle_density = 917
        n_in = n_ice
        k_in = k_ice
        n_ext = n_air
        if small_sizes:
            sz_min = 5e-08
            sz_max = 1e-02
            sz_nbr = 5000
        else:
            sz_min = 2e-05
            sz_max = 1e-01
            sz_nbr = 500

    if water_grain:
        path_to_save_ops = "../Data/OP_data/480band/water_spherical_grains"
        filename = "water_grain"
        medium_type = "air_stp"
        particle_type = "water_grain"
        particle_density = 1000
        n_in = n_water
        k_in = k_water
        n_ext = n_air
        if small_sizes:
            sz_min = 5e-08
            sz_max = 1e-02
            sz_nbr = 5000
        else:
            sz_min = 2e-05
            sz_max = 1e-01
            sz_nbr = 500

    return wvl, ref_index_ice, ref_index_water, n_water, k_water, n_ice, k_ice, n_air, k_in, n_ext, sz_min, sz_max, sz_nbr, path_to_save_temp, air_bbl, water_bbl, ice_grain, water_grain, path_to_save_ops, filename, medium_type, particle_type, particle_density, description

##############################################################################
# Functions
##############################################################################


def net_cdf_out(
    wvl,
    rds,
    particle_density,
    asymmetry,
    ssa,
    ext_xsc,
    sca_xsc,
    abs_xsc,
    ext_cff_vlm,
    sca_cff_vlm,
    abs_cff_vlm,
    ext_cff_mss,
    sca_cff_mss,
    abs_cff_mss,
    eff_radius_analytic,
    eff_radius_resolved,
    median_radius,
    path_to_save,
    descrip,
    medium_tp,
    particle_tp,
    filename,
):

    """
    This function saves the optical properties of a grain/bubble distribution
    to be used in the SNICAR 480bands model.
    The single scattering properties are pre-computed with Mie theory using
    the miepython solver created by Scott Prahl's from the algorithm
    of Wiscombe 1979, and then mixed using a lognormal distribution
    following Flanner et al. 2021.
    """

    file = xr.Dataset(
        data_vars=dict(
            asm_prm=(["wvl"], asymmetry),
            ss_alb=(["wvl"], ssa),
            ext_xsc=(["wvl"], ext_xsc),
            sca_xsc=(["wvl"], sca_xsc),
            abs_xsc=(["wvl"], abs_xsc),
            ext_cff_mss=(["wvl"], ext_cff_mss),
            sca_cff_mss=(["wvl"], sca_cff_mss),
            abs_cff_mss=(["wvl"], abs_cff_mss),
            ext_cff_vlm=(["wvl"], ext_cff_vlm),
            sca_cff_vlm=(["wvl"], sca_cff_vlm),
            abs_cff_vlm=(["wvl"], abs_cff_vlm),
            rds_swa=eff_radius_analytic,
            rds_swr=eff_radius_resolved,
            rds_nma=median_radius,
            gsd=1.5,
            prt_dns=particle_density,
            bnd_nbr=1,
        ),
        coords=dict(wvl=wvl, rds=rds),
        attrs=dict(
            description=descrip,
            medium_type=medium_tp,
            particle_type=particle_tp,
            size_grid="log",
            psd_type="lgn",
        ),
    )

    file["wvl"].attrs = {"long_name": "band center wavelength", "units": "meters"}
    file["rds"].attrs = {"long_name": "particle radius", "units": "meters"}
    file["asm_prm"].attrs = {"long_name": "asymmetry parameter", "units": "fraction"}
    file["ss_alb"].attrs = {
        "long_name": "single scattering albedo",
        "units": "fraction",
    }
    file["ext_xsc"].attrs = {"long_name": "extinction cross-section", "units": "m2"}
    file["sca_xsc"].attrs = {"long_name": "scattering cross-section", "units": "m2"}
    file["abs_xsc"].attrs = {"long_name": "absorption cross-section", "units": "m2"}
    file["ext_cff_mss"].attrs = {
        "long_name": "mass extinction cross-section",
        "units": "m2 kg-1",
    }
    file["sca_cff_mss"].attrs = {
        "long_name": "mass scattering cross-section",
        "units": "m2 kg-1",
    }
    file["abs_cff_mss"].attrs = {
        "long_name": "mass absorption cross-section",
        "units": "m2 kg-1",
    }
    file["ext_cff_vlm"].attrs = {
        "long_name": "volume absorption cross-section",
        "units": "m2 m-3",
    }
    file["sca_cff_vlm"].attrs = {
        "long_name": "volume extinction cross-section",
        "units": "m2 m-3",
    }
    file["abs_cff_vlm"].attrs = {
        "long_name": "volume absorption cross-section",
        "units": "m2 m-3",
    }
    file["rds_swa"].attrs = {
        "long_name": "Analytic surface area-weighted (effective) radius",
        "units": "meters",
    }
    file["rds_swr"].attrs = {
        "long_name": "Resolved surface area-weighted (effective) radius",
        "units": "meters",
    }
    file["rds_nma"].attrs = {
        "long_name": "Analytic number-median radius",
        "units": "meters",
    }
    file["gsd"].attrs = {
        "long_name": "geometric standard deviation of lognormal distribution",
        "units": "unitless",
    }
    file["prt_dns"].attrs = {"long_name": "particle density", "units": "kg m-3"}
    file["bnd_nbr"].attrs = {
        "long_name": "number of bands per wavelength",
        "units": "number",
    }
    file.to_netcdf(
        str(
            f"{path_to_save}/{filename}_{str(int(np.round(eff_radius_analytic*1e6))).zfill(4)}.nc"
        )
    )

    # return file


def n(r, re):
    """

    Parameters
    ----------
    r : log-spaced array of floats (meters)
        radius r
    re : float (meters)
        effective radius

    Returns
    -------
    n : array
        number-distribution
    rn : float
        median radius of the distribution

    """

    # set up lognormal distribution from Flanner et al. 2021
    sigma = 1.5
    rn = re * np.exp(-5 / 2 * np.log(sigma) ** 2)
    coeff = 1 / (np.sqrt(2 * np.pi) * r * np.log(sigma))
    exp = np.exp(-1 / 2 * ((np.log(r) - np.log(rn)) / np.log(sigma)) ** 2)
    n = coeff * exp

    # M. Flanner's code: Apply correction for non-linear spacing in log-grid.
    # This allows us to use N as a pdf, so integrals can be performed
    # as: sum(N*x), rather than sum(N*x*dxm).

    dxm = np.zeros(len(r))
    for i in range(1, len(r) - 1):
        dxm[i] = (r[i + 1] - r[i - 1]) / 2
    dxm[0] = r[1] - r[0]
    dxm[-1] = r[-1] - r[-2]

    return n * dxm, rn


##############################################################################
# Compute OPs of single-size spheres
##############################################################################
def compute_ops_of_single_sized_spheres():
    rds = np.logspace(np.log10(sz_min), np.log10(sz_max), sz_nbr)

    start = time.time()
    for i, radius in enumerate(rds):
        qext, qsca, qback, g = mie.ez_mie(n_in - 1j * k_in, 2 * radius, wvl, n_ext)
        print(f"{np.round(i / len(rds) * 100, 2)} % DONE")

        file = np.concatenate(([qext], [g], [qsca]), axis=0)
        np.save(f"{path_to_save_temp}/{filename}_{radius*1e6}", file)

    end = time.time()
    print(end - start)

    op_filenames = glob.glob(f"{path_to_save_temp}/*.npy")
    radii_ordered = [float(f.split("_")[-1][:-4]) for f in op_filenames]
    op_filenames_ordered = list(
        pd.DataFrame(op_filenames, index=radii_ordered).sort_index().values.flatten()
    )
    op_data = np.array([np.load(fname) for fname in op_filenames_ordered])

    return op_data

##############################################################################
# Compute OPs of lognormal distributions of spheres and save netcdf files
##############################################################################
def compute_ops_of_lognormal_distributions_of_spheres(op_data):
    for re in np.concatenate(
        [np.arange(10e-6, 100e-6, 5e-6), np.arange(100e-6, 5000e-6, 10e-6)]
    ):

        qext = op_data[:, 0, :]
        asm_prm = op_data[:, 1, :]
        qsca = op_data[:, 2, :]

        number_distribution, rn = n(rds, re)

        resolved_re = np.sum(rds**3 * number_distribution) / np.sum(
            rds**2 * number_distribution
        )

        mass_distribution = number_distribution * rds**3
        mass_distribution = mass_distribution / np.sum(mass_distribution)

        qext_wvl = qext.T
        qsca_wvl = qsca.T
        asm_prm_wvl = asm_prm.T

        ext_xsc_wvl = qext_wvl * np.pi * (rds**2)
        sca_xsc_wvl = qsca_wvl * np.pi * (rds**2)

        ext_cff_vlm_wvl = 0.75 * qext_wvl / rds
        sca_cff_vlm_wvl = 0.75 * qsca_wvl / rds

        ext_cff_mss_wvl = 0.75 * qext_wvl / (particle_density * rds)
        sca_cff_mss_wvl = 0.75 * qsca_wvl / (particle_density * rds)

        ext_xsc_out = np.sum(ext_xsc_wvl * number_distribution, axis=1)
        sca_xsc_out = np.sum(sca_xsc_wvl * number_distribution, axis=1)
        abs_xsc_out = ext_xsc_out - sca_xsc_out
        ss_alb_out = sca_xsc_out / ext_xsc_out

        ext_cff_mss_out = np.sum(ext_cff_mss_wvl * mass_distribution, axis=1)
        sca_cff_mss_out = np.sum(sca_cff_mss_wvl * mass_distribution, axis=1)
        abs_cff_mss_out = ext_cff_mss_out - sca_cff_mss_out

        ext_cff_vlm_out = np.sum(ext_cff_vlm_wvl * mass_distribution, axis=1)
        sca_cff_vlm_out = np.sum(sca_cff_vlm_wvl * mass_distribution, axis=1)
        abs_cff_vlm_out = ext_cff_vlm_out - sca_cff_vlm_out

        asm_prm_out = np.sum(
            (asm_prm_wvl * sca_cff_mss_wvl) * mass_distribution, axis=1
        ) / np.sum(sca_cff_mss_wvl * mass_distribution, axis=1)

        net_cdf_out(
            wvl,
            rds,
            particle_density,
            asm_prm_out,
            ss_alb_out,
            ext_xsc_out,
            sca_xsc_out,
            abs_xsc_out,
            ext_cff_vlm_out,
            sca_cff_vlm_out,
            abs_cff_vlm_out,
            ext_cff_mss_out,
            sca_cff_mss_out,
            abs_cff_mss_out,
            re,
            resolved_re,
            rn,
            path_to_save_ops,
            description,
            medium_type,
            particle_type,
            filename,
        )

    return
