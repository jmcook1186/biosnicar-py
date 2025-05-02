#!/usr/bin/python
"""Calculates optical properties of ice/snow column from inputs.

These functions take the user-defined inputs for the physical
properties of the snow and ice and any impurities and calculate
the single scattering albedo, asymmetry parameter and optical
thickness which are then passed to one or other of our radiative
transfer solvers (Toon solver or adding-doubling solver).

"""


import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import pchip

import biosnicar.optical_properties.mie_coated_water_spheres as wcs


def get_layer_OPs(ice, model_config):
    """Calculates optical properties (tauy, ssa, g) of ice column.

    Takes configuration from ice and model_config and uses the data
    to calculate the optical properties of the ice column. There are
    separate routes for layers with granular ice and solid ice.
    Function calls are made to add liquid water coatings or adjust
    the optical properties for aspherical grains where toggled.

    Args:
        ice: instance of Ice class
        model_config: instance of ModelConfig class

    Returns:
        ssa_snw: single scatterign albedo of each layer
        g_snw: asymmetry parameter of each layer
        mac_snw: mass absorption coefficient of each layer

    """

    ssa_snw = np.empty([ice.nbr_lyr, model_config.nbr_wvl])
    mac_snw = np.empty([ice.nbr_lyr, model_config.nbr_wvl])
    g_snw = np.empty([ice.nbr_lyr, model_config.nbr_wvl])
    abs_cff_mss_ice = np.empty(model_config.nbr_wvl)

    # calculations of ice OPs in each layer
    for i in np.arange(0, ice.nbr_lyr, 1):

        if ice.layer_type[i] == 0:  # granular layer

            if ice.shp[i] == 4:  # large hex prisms (geometric optics)
                file_ice = str(
                    model_config.dir_base
                    + model_config.hex_ice_path
                    + ice.op_dir
                    + "{}_{}.nc".format(
                        str(ice.hex_side[i]).rjust(4, "0"), str(ice.hex_length[i])
                    )
                )

            elif ice.shp[i] < 4:
                file_ice = str(
                    model_config.dir_base
                    + model_config.sphere_ice_path
                    + ice.op_dir
                    + "{}.nc".format(str(ice.rds[i]).rjust(4, "0"))
                )

            model_config.file_ice = file_ice

            # if liquid water coatings are applied
            if ice.water[i] > ice.rds[i]:
                ssa_snw[i, :], g_snw[i,:], mac_snw[i,:] = add_water_coating(
                    ice, model_config, ssa_snw[i,:], g_snw[i,:], mac_snw[i,:], i
                )

            else:

                with xr.open_dataset(file_ice) as temp:
                    ssa = temp["ss_alb"].values
                    ssa_snw[i, :] = ssa
                    ext_cff_mss = temp["ext_cff_mss"].values
                    mac_snw[i, :] = ext_cff_mss
                    asm_prm = temp["asm_prm"].values

                    g_snw[i, :] = asm_prm

                    # Correct g for aspherical particles - He et al.(2017)
                    # Applies only when ice.shp!=0
                    # g_snw asymmetry factor parameterization coefficients
                    # (6 bands) from Table 3 & Eqs. 6-7 in He et al. (2017)
                    # assume same values for 4-5 um band, which leads
                    # to very small biases (<3%)

                    if (ice.shp[i] > 0) & (ice.shp[i] < 4):
                        g_snw = correct_for_asphericity(ice, g_snw, ssa_snw, i, model_config)

        # solid ice layer with air/water inclusions
        
        elif (ice.layer_type[i] == 1) or (ice.layer_type[i] == 2):

            if ice.cdom[i]:
                cdom = pd.read_csv(
                    model_config.dir_base + "Data/OP_data/k_cdom_240_750.csv"
                )
                cdom_ref_idx_im = np.array(cdom).flatten()
                # rescale to SNICAR resolution
                cdom_ref_idx_im_rescaled = cdom_ref_idx_im[::10]
                ice.ref_idx_im[3:54] = np.fmax(
                    ice.ref_idx_im[3:54], cdom_ref_idx_im_rescaled
                )
                
            # neglecting air mass:
            vlm_frac_ice = (ice.rho[i] - ice.lwc[i] * 1000) / 917
            vlm_frac_air = 1 - ice.lwc[i] - vlm_frac_ice

            # get effective radius
            rd = f"{ice.rds[i]}"
            rd = rd.rjust(4, "0")
            file_ice_path = str(model_config.dir_base +
                                model_config.bubbly_ice_path 
                                + "bbl_{}.nc").format(rd)
            file_ice = xr.open_dataset(file_ice_path)
            
            # air bbl ssps
            sca_cff_vlm_air_bbl = file_ice["sca_cff_vlm"].values
            g_air_bbl = file_ice["asm_prm"].values

            if ice.lwc[i] == 0:
            
                abs_cff_mss_ice[:] = ((4 * np.pi * ice.ref_idx_im) / (model_config.wavelengths * 1e-6)) / 917
                mac_snw[i, :] = (
                    (sca_cff_vlm_air_bbl * vlm_frac_air) / ice.rho[i]
                ) + abs_cff_mss_ice

                ssa_snw[i, :] = (
                    (sca_cff_vlm_air_bbl * vlm_frac_air) / ice.rho[i]
                ) / mac_snw[i, :]
                
                g_snw[i, :] = g_air_bbl
            
            elif ice.lwc[i] != 0:

                # water bubbles ssps
                file_water = xr.open_dataset(
                    str(model_config.dir_base + model_config.bubbly_ice_path + "bbl_water_{}.nc").format(rd)
                )
                sca_cff_vlm_water = file_water["sca_cff_vlm"].values
                ext_cff_vlm_water = file_water["ext_cff_vlm"].values
                g_water = file_water["asm_prm"].values

                vlm_frac_lw_in_ice = ice.lwc[i] * (1 - ice.lwc_pct_bbl)
                vlm_frac_lw_in_bbl = ice.lwc[i] * ice.lwc_pct_bbl

                # neglecting air absorption:
                abs_cff_mss_ice[:] = (vlm_frac_ice * 917 / ice.rho[i]) * (
                    4 * np.pi * ice.ref_idx_im / (model_config.wavelengths * 1e-6)
                ) / 917 + (vlm_frac_lw_in_ice * 1000 / ice.rho[i]) * (
                    4 * np.pi * ice.ref_idx_im_water / (model_config.wavelengths * 1e-6)
                ) / 1000

                # volume weighted assymetry parameter
                g_snw[i, :] = (g_air_bbl * vlm_frac_air + g_water * vlm_frac_lw_in_bbl) / (
                    vlm_frac_lw_in_bbl + vlm_frac_air
                )

                # volume weighted extinction coefficient
                mac_snw[i, :] = (
                    (sca_cff_vlm_air_bbl * vlm_frac_air) / ice.rho[i]
                    + (ext_cff_vlm_water * vlm_frac_lw_in_bbl) / ice.rho[i]
                ) + abs_cff_mss_ice

                # volume weighted scattering coefficient
                ssa_snw[i, :] = (
                    (sca_cff_vlm_air_bbl * vlm_frac_air) / ice.rho[i]
                    + (sca_cff_vlm_water * vlm_frac_lw_in_bbl) / ice.rho[i]
                ) / mac_snw[i, :]
                
                
        # granular layer with mixed ice and water spheres
        elif ice.layer_type[i] == 3:

            # neglecting air mass:
            vlm_frac_ice = (ice.rho[i] - ice.lwc[i] * 1000) / 917
            vlm_frac_air = 1 - ice.lwc[i] - vlm_frac_ice

            # get effective radius
            rd = f"{ice.rds[i]}"
            rd = rd.rjust(4, "0")
            ssps_ice = xr.open_dataset(
                str(model_config.dir_base +
                    model_config.sphere_ice_path
                    + ice.op_dir
                    + "{}.nc".format(str(ice.rds[i]).rjust(4, "0"))
                )
            )
            ssps_water = xr.open_dataset(
                str(model_config.dir_base +
                    model_config.sphere_water_path
                    + "water_grain_{}.nc".format(str(ice.rds[i]).rjust(4, "0"))
                )
            )

            g = (
                ssps_water["asm_prm"].values * ice.lwc[i]
                + ssps_ice["asm_prm"].values * vlm_frac_ice
            ) / (ice.lwc[i] + vlm_frac_ice)

            g_snw[i, :] = g

            ext_cff_mss = (
                ssps_water["ext_cff_vlm"].values * ice.lwc[i]
                + ssps_ice["ext_cff_vlm"].values * vlm_frac_ice
            ) / (ice.rho[i])

            mac_snw[i, :] = ext_cff_mss

            ssa = (
                (
                    ssps_water["sca_cff_vlm"].values * ice.lwc[i]
                    + ssps_ice["sca_cff_vlm"].values * vlm_frac_ice
                )
                / ice.rho[i]
                / mac_snw[i, :]
            )

            ssa_snw[i, :] = ssa


    return ssa_snw, g_snw, mac_snw


def add_water_coating(ice, model_config, ssa_snw, g_snw, mac_snw, i):

    """Recalculates layer optical properties where grains are coated in liquid water.

    Feature originally added by Niklas Bohn. Where value of water exceeds value of rds
    for a given layer it is interpreted as having a liquid water film. in this case
    Mie calculations for a coated sphere are executed with the outer coating havign radius
    water - rds.

    Args:
        ice: instance of Ice class
        model_config: instance of ModelConfig class
        ssa_snw: single scattering albedo of each layer
        g_snw: asymmetry parameter for each layer
        mac_snw: mass absorption coefficient of each layer
        i: layer counter

    Returns:
        ssa_snw: updated single scattering albedo for each layer
        g_snw: updated asymmetry parameter for each layer
        mac_snw: updated mass absorption coefficient for each layer

    Raises:
        ValueError if ice.shp!= 0 (i.e. grains not spherical)

    """

    if ice.shp[i] != 0:
        raise ValueError("Water coating can only be applied to spheres")

    res = wcs.miecoated_driver(
        rice=ice.rds[i],
        rwater=ice.water[i],
        fn_ice=model_config.dir_base+model_config.fn_ice,
        rf_ice=ice.rf,
        fn_water=model_config.dir_base+model_config.fn_water,
        wvl=model_config.wavelengths,
    )


    ssa_snw = res["ssa"]
    g_snw = res["asymmetry"]

    with xr.open_dataset(model_config.file_ice) as temp:

        ext_cff_mss = temp["ext_cff_mss"].values
        mac_snw = ext_cff_mss

    return ssa_snw, g_snw, mac_snw


def correct_for_asphericity(ice, g_snw, ssa_snw, i, model_config):
    """Adjusts asymmetry parameter for aspherical grains.

    Implements work from Fu et al. 2007 and He et al. 2017.
    Asymmetry parameter is adjusted to account for asphericity
    for the defined shape of each layer.

    Ice grain shape can be
    0 = sphere,
    1 = spheroid,
    2 = hexagonal plate,
    3 = koch snowflake,
    4 = hexagonal prisms

    Args:
        ice: instance of Ice class
        g_snw: asymmetry parameter for each layer
        i: layer counter

    Returns:
        g_snw: updated asymmetry parameter for layer
    """

    g_wvl = np.array([0.25, 0.70, 1.41, 1.90, 2.50, 3.50, 4.00, 5.00])

    g_wvl_center = np.array(g_wvl[1:8]) / 2 + np.array(g_wvl[0:7]) / 2
    g_b0 = np.array(
        [
            9.76029e-01,
            9.67798e-01,
            1.00111e00,
            1.00224e00,
            9.64295e-01,
            9.97475e-01,
            9.97475e-01,
        ]
    )

    g_b1 = np.array(
        [
            5.21042e-01,
            4.96181e-01,
            1.83711e-01,
            1.37082e-01,
            5.50598e-02,
            8.48743e-02,
            8.48743e-02,
        ]
    )

    g_b2 = np.array(
        [
            -2.66792e-04,
            1.14088e-03,
            2.37011e-04,
            -2.35905e-04,
            8.40449e-04,
            -4.71484e-04,
            -4.71484e-04,
        ]
    )

    # Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007
    g_f07_c2 = np.array(
        [
            1.349959e-1,
            1.115697e-1,
            9.853958e-2,
            5.557793e-2,
            -1.233493e-1,
            0.0,
            0.0,
        ]
    )
    g_f07_c1 = np.array(
        [
            -3.987320e-1,
            -3.723287e-1,
            -3.924784e-1,
            -3.259404e-1,
            4.429054e-2,
            -1.726586e-1,
            -1.726586e-1,
        ]
    )
    g_f07_c0 = np.array(
        [
            7.938904e-1,
            8.030084e-1,
            8.513932e-1,
            8.692241e-1,
            7.085850e-1,
            6.412701e-1,
            6.412701e-1,
        ]
    )
    g_f07_p2 = np.array(
        [
            3.165543e-3,
            2.014810e-3,
            1.780838e-3,
            6.987734e-4,
            -1.882932e-2,
            -2.277872e-2,
            -2.277872e-2,
        ]
    )
    g_f07_p1 = np.array(
        [
            1.140557e-1,
            1.143152e-1,
            1.143814e-1,
            1.071238e-1,
            1.353873e-1,
            1.914431e-1,
            1.914431e-1,
        ]
    )
    g_f07_p0 = np.array(
        [
            5.292852e-1,
            5.425909e-1,
            5.601598e-1,
            6.023407e-1,
            6.473899e-1,
            4.634944e-1,
            4.634944e-1,
        ]
    )

    fs_hex = 0.788  # shape factor for hex plate

    # eff grain diameter
    diam_ice = 2.0 * ice.rds[i] / 0.544

    if ice.shp_fctr[i] == 0:
        # default shape factor for koch snowflake;
        # He et al. (2017), Table 1
        fs_koch = 0.712

    else:

        fs_koch = ice.shp_fctr[i]

    if ice.grain_ar[i] == 0:
        # default aspect ratio for koch
        # snowflake; He et al. (2017), Table 1
        ar_tmp = 2.5

    else:

        ar_tmp = ice.grain_ar[i]

    # Eq.7, He et al. (2017)
    g_snw_cg_tmp = g_b0 * (fs_koch / fs_hex) ** g_b1 * diam_ice**g_b2

    # Eqn. 3.3 in Fu (2007)
    gg_snw_f07_tmp = (
        g_f07_p0 + g_f07_p1 * np.log(ar_tmp) + g_f07_p2 * (np.log(ar_tmp)) ** 2
    )

    # 1 = spheroid, He et al. (2017)
    if ice.shp[i] == 1:

        # effective snow grain diameter
        diam_ice = 2.0 * ice.rds[i]

        # default shape factor for spheroid;
        # He et al. (2017), Table 1
        if ice.shp_fctr[i] == 0:

            fs_sphd = 0.929

        else:
            # if shp_factor not 0,
            # then use user-defined value
            fs_sphd = ice.shp_fctr[i]

        if ice.grain_ar[i] == 0:
            # default aspect ratio for spheroid;
            # He et al. (2017), Table 1
            ar_tmp = 0.5

        else:

            ar_tmp = ice.grain_ar[i]

        # Eq.7, He et al. (2017)
        g_snw_cg_tmp = g_b0 * (fs_sphd / fs_hex) ** g_b1 * diam_ice**g_b2

        # Eqn. 3.1 in Fu (2007)
        gg_snw_F07_tmp = g_f07_c0 + g_f07_c1 * ar_tmp + g_f07_c2 * ar_tmp**2

    # 3=hexagonal plate,
    # He et al. 2017 parameterization
    if ice.shp[i] == 2:

        # effective snow grain diameter
        diam_ice = 2.0 * ice.rds[i]

        if ice.shp_fctr[i] == 0:
            # default shape factor for
            # hexagonal plates;
            # He et al. (2017), Table 1
            fs_hex0 = 0.788

        else:

            fs_hex0 = ice.shp_fctr[i]

        if ice.grain_ar[i] == 0:
            # default aspect ratio
            # for hexagonal plate;
            # He et al. (2017), Table 1
            ar_tmp = 2.5

        else:

            ar_tmp = ice.grain_ar[i]

        # Eq.7, He et al. (2017)
        g_snw_cg_tmp = g_b0 * (fs_hex0 / fs_hex) ** g_b1 * diam_ice**g_b2

        # Eqn. 3.3 in Fu (2007)
        gg_snw_F07_tmp = (
            g_f07_p0 + g_f07_p1 * np.log(ar_tmp) + g_f07_p2 * (np.log(ar_tmp)) ** 2
        )

    # 4=koch snowflake,
    # He et al. (2017)
    #  parameterization
    if ice.shp[i] == 3:

        # effective snow grain diameter
        diam_ice = 2.0 * ice.rds[i] / 0.544

        if ice.shp_fctr[i] == 0:
            # default shape factor
            # for koch snowflake;
            # He et al. (2017), Table 1
            fs_koch = 0.712

        else:

            fs_koch = ice.shp_fctr[i]

        # default aspect ratio for
        # koch snowflake; He et al. (2017), Table 1
        if ice.grain_ar[i] == 0:

            ar_tmp = 2.5

        else:

            ar_tmp = ice.grain_ar[i]

        # Eq.7, He et al. (2017)
        g_snw_cg_tmp = g_b0 * (fs_koch / fs_hex) ** g_b1 * diam_ice**g_b2

        # Eqn. 3.3 in Fu (2007)
        gg_snw_F07_tmp = (
            g_f07_p0 + g_f07_p1 * np.log(ar_tmp) + g_f07_p2 * (np.log(ar_tmp)) ** 2
        )

    # 6 wavelength bands for g_snw to be
    # interpolated into 480-bands of SNICAR
    # shape-preserving piecewise interpolation
    # into 480-bands
    g_Cg_intp = pchip(g_wvl_center, g_snw_cg_tmp)(model_config.wavelengths)
    gg_f07_intp = pchip(g_wvl_center, gg_snw_F07_tmp)(model_config.wavelengths)
    g_snw_F07 = (
        gg_f07_intp + (1.0 - gg_f07_intp) / ssa_snw[i, :] / 2
    )  # Eq.2.2 in Fu (2007)
    # Eq.6, He et al. (2017)
    g_snw[i, :] = g_snw_F07 * g_Cg_intp
    g_snw[i, 381:480] = g_snw[i, 380]
    # assume same values for 4-5 um band,
    # with v small biases (<3%)

    g_snw[g_snw <= 0] = 0.00001
    g_snw[g_snw > 0.99] = 0.99  # avoid unreasonable
    # values (so far only occur in large-size spheroid cases)

    return g_snw


def mix_in_impurities(ssa_snw, g_snw, mac_snw, ice, impurities, model_config):
    """Updates optical properties for the presence of light absorbing particles.

    Takes the optical properties of the clean ice column and adjusts them for
    the presence of light absorbing particles in the ice. Each impurity is an
    instance of the Impurity class whose attributes include the path to the
    specific optical properties for that impurity. Its concentration is generally
    provided in ppb, but concentration of algae can also be given in cells/mL.


    Args:
        ssa_snw: single scattering albedo of eahc layer
        g_snw: asymmetry parameter for each layer
        mac_snw: mass absorption coefficient of each layer
        ice: instance of Ice class
        impurities: array containing instances of Impurity class
        model_config: instance of ModelConfig class

    Returns:
        tau: updated optical thickness
        ssa: updated single scattering albedo
        g: updated asymmetry parameter
        L_snw: mass of ice ine ach layer

    """

    ssa_aer = np.zeros([len(impurities), model_config.nbr_wvl])
    mac_aer = np.zeros([len(impurities), model_config.nbr_wvl])
    g_aer = np.zeros([len(impurities), model_config.nbr_wvl])
    mss_aer = np.zeros([ice.nbr_lyr, len(impurities)])
    g_sum = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    ssa_sum = np.zeros([ice.nbr_lyr, len(impurities), model_config.nbr_wvl])
    tau = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    ssa = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    g = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    L_aer = np.zeros([ice.nbr_lyr, len(impurities)])
    tau_aer = np.zeros([ice.nbr_lyr, len(impurities), model_config.nbr_wvl])
    tau_sum = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    ssa_sum = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    L_snw = np.zeros(ice.nbr_lyr)
    tau_snw = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])

    for i, impurity in enumerate(impurities):

        g_aer[i, :] = impurity.g
        ssa_aer[i, :] = impurity.ssa

        if impurity.unit == 1:

            mss_aer[0 : ice.nbr_lyr, i] = (
                np.array(impurity.conc) / 917 * 10**6
            ) 

        else:
            mss_aer[0 : ice.nbr_lyr, i] = (
                np.array(impurity.conc) * 1e-9
            ) 

    # for each layer, the layer mass (L) is density * layer thickness
    # for each layer the optical ice.depth is
    # the layer mass * the mass extinction coefficient
    # first for the ice in each layer

    for i in range(ice.nbr_lyr):

        L_snw[i] = ice.rho[i] * ice.dz[i]

        for j, impurity in enumerate(impurities):

            mac_aer[j, :] = impurity.mac

            # kg ice m-2 * cells kg-1 ice = cells m-2
            L_aer[i, j] = L_snw[i] * mss_aer[i, j]
            # cells m-2 * m2 cells-1

            tau_aer[i, j, :] = L_aer[i, j] * mac_aer[j, :]
            tau_sum[i, :] = tau_sum[i, :] + tau_aer[i, j, :]
            ssa_sum[i, :] = ssa_sum[i, :] + (tau_aer[i, j, :] * ssa_aer[j, :])
            g_sum[i, :] = g_sum[i, :] + (tau_aer[i, j, :] * ssa_aer[j, :] * g_aer[j, :])

            # ice mass = snow mass - impurity mass (generally tiny correction)
            # if aer == algae and L_aer is in cells m-2, should be converted
            # to m-2 kg-1 : 1 cell = 1ng = 10**(-12) kg

            if impurity.unit == 1:

                L_snw[i] = L_snw[i] - L_aer[i, j] * 10 ** (-12)

            else:
                L_snw[i] = L_snw[i] - L_aer[i, j]

        tau_snw[i, :] = L_snw[i] * mac_snw[i, :]

        # finally, for each layer calculate the effective ssa, tau and g
        # for the snow+LAP
        tau[i, :] = tau_sum[i, :] + tau_snw[i, :]
        ssa[i, :] = (1 / tau[i, :]) * (ssa_sum[i, :] + (ssa_snw[i, :] * tau_snw[i, :]))
        g[i, :] = (1 / (tau[i, :] * (ssa[i, :]))) * (
            g_sum[i, :] + (g_snw[i, :] * ssa_snw[i, :] * tau_snw[i, :])
        )

    # just in case any unrealistic values arise (none detected so far)
    ssa[ssa <= 0] = 0.00000001
    ssa[ssa >= 1] = 0.99999999
    g[g <= 0] = 0.00001
    g[g > 0.99] = 0.99

    return tau, ssa, g, L_snw


if __name__ == "__main__":
    pass
