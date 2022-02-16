import collections as c
import math
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import pchip
import mie_coated_water_spheres as wcs
from classes import *

def get_layer_OPs(ice, impurities, model_cfg, ice_cfg, rtm_cfg):

    ssa_snw = np.empty([ice.nbr_lyr, rtm_cfg["nbr_wvl"]])
    mac_snw = np.empty([ice.nbr_lyr, rtm_cfg["nbr_wvl"]])
    g_snw = np.empty([ice.nbr_lyr, rtm_cfg["nbr_wvl"]])
    abs_cff_mss_ice = np.empty(rtm_cfg["nbr_wvl"])

    # calculations of ice OPs in each layer
    for i in np.arange(0, ice.nbr_lyr, 1):

        if ice.layer_type[i] == 0:  # granular layer

            if ice.shp[i] == 4:  # large hex prisms (geometric optics)
                file_ice = str(
                    model_cfg["PATHS"]["DIR_BASE"] +
                    ice_cfg["PATHS"]["HEX_ICE"]
                    + ice.op_dir
                    + "{}_{}.nc".format(
                        str(ice.hex_side[i]).rjust(4, "0"), str(ice.hex_length[i])
                    )
                )

            elif ice.shp[i] < 4:
                file_ice = str(
                    model_cfg["PATHS"]["DIR_BASE"] +
                    ice_cfg["PATHS"]["SPHERE_ICE"]
                    + ice.op_dir
                    + "{}.nc".format(str(ice.rds[i]).rjust(4, "0"))
                )


            # if liquid water coatings are applied
            if ice.water[i] > ice.rds[i]:
                ssa_snw[i, ice.nbr_wvl], g_snw[i], mac_snw[i] = add_water_coating(ice, ice_cfg, model_cfg, ssa_snw, g_snw, mac_snw, i)

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
                        g_snw = correct_for_asphericity(ice, ice_cfg, model_cfg, g_snw, i)


        else:  # solid ice layer (ice.layer_type == 1)

            if ice.cdom_layer[i]:
                cdom_refidx_im = np.array(
                    pd.read_csv(ice_cfg["PATHS"]["RI_ICE"] + "k_cdom_240_750.csv")
                ).flatten()

                # rescale to SNICAR resolution
                cdom_refidx_im_rescaled = cdom_refidx_im[::10]
                refidx_im[3:54] = np.fmax(refidx_im[3:54], cdom_refidx_im_rescaled)

            rd = f"{ice.rds[i]}".rjust(4, "0")
            file_ice = str(ice_cfg["PATHS"]["BUBBLY_ICE"] + "bbl_{}.nc").format(rd)
            file = xr.open_dataset(file_ice)
            sca_cff_vlm = file["sca_cff_vlm"].values
            g_snw[i, :] = file["asm_prm"].values
            abs_cff_mss_ice[:] = ((4 * np.pi * refidx_im) / (wvl * 1e-6)) / 917
            vlm_frac_air = (917 - ice.rho_layers[i]) / 917
            mac_snw[i, :] = (
                (sca_cff_vlm * vlm_frac_air) / ice.rho_layers[i]
            ) + abs_cff_mss_ice
            ssa_snw[i, :] = ((sca_cff_vlm * vlm_frac_air) / ice.rho_layers[i]) / mac_snw[
                i, :
            ]

    return ssa_snw, g_snw, mac_snw



def add_water_coating(ice, ice_cfg, model_cfg, ssa_snw, g_snw, mac_snw, i):

    """
    adds liquid water coating to spherical ice grains

    """

    if ice.shp[i] != 0:
        raise ValueError("Water coating can only be applied to spheres")

    fn_ice = model_cfg["PATHS"]["DIR_BASE"] + ice_cfg["PATHS"]["RFX_ICE"]
    fn_water = model_cfg["PATHS"]["DIR_BASE"] + ice_cfg["PATHS"]["RFX_WATER"]

    res = wcs.miecoated_driver(
        rice=ice.rds[i],
        water=ice.water[i],
        fn_ice=fn_ice,
        rf_ice=ice.rf_ice,
        fn_water=fn_water,
        wvl = get_wavelengths(),
    )

    ssa_snw[i, :] = res["ssa"]
    g_snw[i, :] = res["asymmetry"]

    with xr.open_dataset(file_ice) as temp:
        ext_cff_mss = temp["ext_cff_mss"].values
        mac_snw[i, :] = ext_cff_mss
    
    return ssa_snw, g_snw, mac_snw


def correct_for_asphericity(ice, ice_cfg, model_cfg, g_snw, i):
    
    g_wvl = np.array(
    [0.25, 0.70, 1.41, 1.90, 2.50, 3.50, 4.00, 5.00])
    
    g_wvl_center = (np.array(g_wvl[1:8]) / 2 + np.array(g_wvl[0:7]) / 2)
    g_b0 = np.array([
        9.76029e-01,
        9.67798e-01,
        1.00111e00,
        1.00224e00,
        9.64295e-01,
        9.97475e-01,
        9.97475e-01]
    )
    
    g_b1 = np.array([
        5.21042e-01,
        4.96181e-01,
        1.83711e-01,
        1.37082e-01,
        5.50598e-02,
        8.48743e-02,
        8.48743e-02]
    )
    
    g_b2 = np.array([
        -2.66792e-04,
        1.14088e-03,
        2.37011e-04,
        -2.35905e-04,
        8.40449e-04,
        -4.71484e-04,
        -4.71484e-04,]
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
    g_snw_cg_tmp = (
        g_b0 * (fs_koch / fs_hex) ** g_b1 * diam_ice ** g_b2
    )

    # Eqn. 3.3 in Fu (2007)
    gg_snw_f07_tmp = (
        g_f07_p0
        + g_f07_p1 * np.log(ar_tmp)
        + g_f07_p2 * (np.log(ar_tmp)) ** 2
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
        g_snw_cg_tmp = (
            g_b0 * (fs_sphd / fs_hex) ** g_b1 * diam_ice ** g_b2
        )

        # Eqn. 3.1 in Fu (2007)
        gg_snw_F07_tmp = (
            g_f07_c0 + g_f07_c1 * ar_tmp + g_f07_c2 * ar_tmp ** 2
        )

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
        g_snw_cg_tmp = (
            g_b0 * (fs_hex0 / fs_hex) ** g_b1 * diam_ice ** g_b2
        )

        # Eqn. 3.3 in Fu (2007)
        gg_snw_F07_tmp = (
            g_f07_p0
            + g_f07_p1 * np.log(ar_tmp)
            + g_f07_p2 * (np.log(ar_tmp)) ** 2
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
        g_snw_cg_tmp = (
            g_b0 * (fs_koch / fs_hex) ** g_b1 * diam_ice ** g_b2
        )

        # Eqn. 3.3 in Fu (2007)
        gg_snw_F07_tmp = (
            g_f07_p0
            + g_f07_p1 * np.log(ar_tmp)
            + g_f07_p2 * (np.log(ar_tmp)) ** 2
        )

    # 6 wavelength bands for g_snw to be
    # interpolated into 480-bands of SNICAR
    # shape-preserving piecewise interpolation
    # into 480-bands
    g_Cg_intp = pchip(g_wvl_center, g_snw_cg_tmp)(wvl)
    gg_f07_intp = pchip(g_wvl_center, gg_snw_F07_tmp)(wvl)
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