import sys

sys.path.append("./src")
import collections as c
import numpy as np
from snicar_feeder import snicar_feeder


def call_snicar(params):
    """
    enables snicar to be called as a function rather than by
    running the driver script
    """

    inputs = c.namedtuple(
        "inputs",
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

    ##############################
    ## 2) Set working directory
    ##############################

    # set dir_base to the location of the BioSNICAR_GO_PY folder
    inputs.dir_base = "/home/joe/Code/BioSNICAR_GO_PY/"
    inputs.verbosity = 0  # 1 to print real-time updates

    #######################################
    ## 4) RADIATIVE TRANSFER CONFIGURATION
    #######################################

    inputs.direct = 1
    inputs.aprx_typ = 1
    inputs.delta = 1
    inputs.solzen = params.solzen

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
    inputs.incoming_i = 4

    ###############################################################
    ## 4) SET UP ICE/SNOW LAYERS
    # For granular layers only, choose toon
    # For granular layers + Fresnel layers below, choose add_double
    ###############################################################

    inputs.toon = False
    inputs.add_double = True

    inputs.dz = params.dz
    inputs.nbr_lyr = len(params.dz)
    inputs.layer_type = params.layer_type
    inputs.cdom_layer = [0, 0]
    inputs.rho_layers = params.rho_layers
    inputs.nbr_wvl = 480
    inputs.R_sfc = np.genfromtxt(
        "./Data/OP_data/480band/rain_polished_ice_spectrum.csv", delimiter="csv"
    )  # import underlying ice from file

    inputs.rf_ice = 2
    inputs.grain_shp = [0] * len(params.dz)
    inputs.grain_rds = params.grain_rds
    inputs.rwater = [0] * len(params.dz)
    inputs.side_length = 0
    inputs.depth = 0
    inputs.shp_fctr = [0] * len(params.dz)
    inputs.grain_ar = [0] * len(params.dz)

    #######################################
    ## 5) SET LAP CHARACTERISTICS
    #######################################

    # Define total number of different LAPs/aerosols in model
    inputs.nbr_aer = 30

    # define units for algae absorption cross section input file
    # 0 = m2/kg for MAC, ppb for mss_cnc (this is default)
    # 1 = m2/cell for MAC, cell/mL for mss_cnc
    inputs.GA_units = 1  # glacier algae
    inputs.SA_units = 1  # snow algae

    # determine C_factor (can be None or a number)
    # this is the concentrating factor that accounts for
    # resolution difference in field samples and model layers
    inputs.c_factor_GA = params.c_factor_GA
    inputs.c_factor_SA = 0

    # Set names of files containing the optical properties of these LAPs:
    inputs.file_soot1 = "mie_sot_ChC90_dns_1317.nc"
    inputs.file_soot2 = "miecot_slfsot_ChC90_dns_1317.nc"
    inputs.file_brwnC1 = "brC_Kirch_BCsd.nc"
    inputs.file_brwnC2 = "brC_Kirch_BCsd_slfcot.nc"
    inputs.file_dust1 = "dust_balkanski_central_size1.nc"
    inputs.file_dust2 = "dust_balkanski_central_size2.nc"
    inputs.file_dust3 = "dust_balkanski_central_size3.nc"
    inputs.file_dust4 = "dust_balkanski_central_size4.nc"
    inputs.file_dust5 = "dust_balkanski_central_size5.nc"
    inputs.file_ash1 = "volc_ash_eyja_central_size1.nc"
    inputs.file_ash2 = "volc_ash_eyja_central_size2.nc"
    inputs.file_ash3 = "volc_ash_eyja_central_size3.nc"
    inputs.file_ash4 = "volc_ash_eyja_central_size4.nc"
    inputs.file_ash5 = "volc_ash_eyja_central_size5.nc"
    inputs.file_ash_st_helens = "volc_ash_mtsthelens_20081011.nc"
    inputs.file_Skiles_dust1 = "dust_skiles_size1.nc"
    inputs.file_Skiles_dust2 = "dust_skiles_size2.nc"
    inputs.file_Skiles_dust3 = "dust_skiles_size3.nc"
    inputs.file_Skiles_dust4 = "dust_skiles_size4.nc"
    inputs.file_Skiles_dust5 = "dust_skiles_size5.nc"
    inputs.file_GreenlandCentral1 = "dust_greenland_central_size1.nc"
    inputs.file_GreenlandCentral2 = "dust_greenland_central_size2.nc"
    inputs.file_GreenlandCentral3 = "dust_greenland_central_size3.nc"
    inputs.file_GreenlandCentral4 = "dust_greenland_central_size4.nc"
    inputs.file_GreenlandCentral5 = "dust_greenland_central_size5.nc"
    inputs.file_Cook_Greenland_dust_L = "dust_greenland_Cook_LOW_20190911.nc"
    inputs.file_Cook_Greenland_dust_C = "dust_greenland_Cook_CENTRAL_20190911.nc"
    inputs.file_Cook_Greenland_dust_H = "dust_greenland_Cook_HIGH_20190911.nc"
    inputs.file_snw_alg = "snw_alg_r025um_chla020_chlb025_cara150_carb140.nc"
    inputs.file_glacier_algae = "GA_Chevrollier2022_r4.9_L18.8.nc"

    inputs.mss_cnc_soot1 = [0] * len(params.dz)
    inputs.mss_cnc_soot2 = [0] * len(params.dz)
    inputs.mss_cnc_brwnC1 = [0] * len(params.dz)
    inputs.mss_cnc_brwnC2 = [0] * len(params.dz)
    inputs.mss_cnc_dust1 = [0] * len(params.dz)
    inputs.mss_cnc_dust2 = [0] * len(params.dz)
    inputs.mss_cnc_dust3 = [0] * len(params.dz)
    inputs.mss_cnc_dust4 = [0] * len(params.dz)
    inputs.mss_cnc_dust5 = [0] * len(params.dz)
    inputs.mss_cnc_ash1 = [0] * len(params.dz)
    inputs.mss_cnc_ash2 = [0] * len(params.dz)
    inputs.mss_cnc_ash3 = [0] * len(params.dz)
    inputs.mss_cnc_ash4 = [0] * len(params.dz)
    inputs.mss_cnc_ash5 = [0] * len(params.dz)
    inputs.mss_cnc_ash_st_helens = [0] * len(params.dz)
    inputs.mss_cnc_Skiles_dust1 = [0] * len(params.dz)
    inputs.mss_cnc_Skiles_dust2 = [0] * len(params.dz)
    inputs.mss_cnc_Skiles_dust3 = [0] * len(params.dz)
    inputs.mss_cnc_Skiles_dust4 = [0] * len(params.dz)
    inputs.mss_cnc_Skiles_dust5 = [0] * len(params.dz)
    inputs.mss_cnc_GreenlandCentral1 = [0] * len(params.dz)
    inputs.mss_cnc_GreenlandCentral2 = [0] * len(params.dz)
    inputs.mss_cnc_GreenlandCentral3 = [0] * len(params.dz)
    inputs.mss_cnc_GreenlandCentral4 = [0] * len(params.dz)
    inputs.mss_cnc_GreenlandCentral5 = [0] * len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_L = [0] * len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_C = [0] * len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_H = [0] * len(params.dz)
    inputs.mss_cnc_snw_alg = [0] * len(params.dz)
    inputs.mss_cnc_glacier_algae = params.mss_cnc_glacier_algae

    outputs = snicar_feeder(inputs)

    return outputs.albedo, outputs.BBA
