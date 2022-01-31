import collections

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from snicar_feeder import snicar_feeder


def generate_snicar_params(
    layer_type, density, dz, alg, solzen, reff, bc, dust1, dust5
):

    rho_layers = [density] * len(dz)
    reff = [reff] * len(dz)
    layer_type = [layer_type] * len(dz)
    mss_cnc_glacier_algae = [alg] * len(dz)
    mss_cnc_soot1 = [0] * len(dz)
    mss_cnc_soot1[0] = bc
    mss_cnc_dust1 = [0] * len(dz)
    mss_cnc_dust1[0] = dust1
    mss_cnc_dust5 = [0] * len(dz)
    mss_cnc_dust5[0] = dust5

    params = collections.namedtuple(
        "params",
        "rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, mss_cnc_soot1, mss_cnc_dust1, mss_cnc_dust5, solzen",
    )
    params.grain_rds = reff
    params.rho_layers = rho_layers
    params.layer_type = layer_type
    params.dz = dz
    params.mss_cnc_glacier_algae = mss_cnc_glacier_algae
    params.solzen = solzen
    params.mss_cnc_soot1 = mss_cnc_soot1
    params.mss_cnc_dust1 = mss_cnc_dust1
    params.mss_cnc_dust5 = mss_cnc_dust5

    return params


def call_snicar(params):

    inputs = collections.namedtuple(
        "inputs",
        [
            "dir_base",
            "verbosity",
            "rf_ice",
            "incoming_i",
            "DIRECT",
            "layer_type",
            "APRX_TYP",
            "DELTA",
            "solzen",
            "TOON",
            "ADD_DOUBLE",
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
            "SA_units",
            "GA_units",
            "Cfactor_SA",
            "Cfactor_GA",
            "cdom_layer",
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
            "FILE_soot1",
            "FILE_soot2",
            "FILE_brwnC1",
            "FILE_brwnC2",
            "FILE_dust1",
            "FILE_dust2",
            "FILE_dust3",
            "FILE_dust4",
            "FILE_dust5",
            "FILE_ash1",
            "FILE_ash2",
            "FILE_ash3",
            "FILE_ash4",
            "FILE_ash5",
            "FILE_ash_st_helens",
            "FILE_Skiles_dust1",
            "FILE_Skiles_dust2",
            "FILE_Skiles_dust3",
            "FILE_Skiles_dust4",
            "FILE_Skiles_dust5",
            "FILE_GreenlandCentral1",
            "FILE_GreenlandCentral2",
            "FILE_GreenlandCentral3",
            "FILE_GreenlandCentral4",
            "FILE_GreenlandCentral5",
            "FILE_Cook_Greenland_dust_L",
            "FILE_Cook_Greenland_dust_C",
            "FILE_Cook_Greenland_dust_H",
            "FILE_snw_alg",
            "FILE_glacier_algae",
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
    inputs.verbosity = 0
    savepath = inputs.dir_base  # base path for saving figures

    ################################
    ## 3) Choose plot/print options
    ################################

    show_figs = False  # toggle to display spectral albedo figure
    save_figs = False  # toggle to save spectral albedo figure to file
    print_BBA = False  # toggle to print broadband albedo to terminal
    print_band_ratios = False  # toggle to print various band ratios to terminal
    smooth = False  # apply optional smoothing function (Savitzky-Golay filter)
    window_size = 9  # if applying smoothing filter, define window size
    poly_order = 3  # if applying smoothing filter, define order of polynomial

    #######################################
    ## 4) RADIATIVE TRANSFER CONFIGURATION
    #######################################

    inputs.DIRECT = 1  # 1= Direct-beam incident flux, 0= Diffuse incident flux
    inputs.APRX_TYP = 1  # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    inputs.DELTA = 1  # 1= Apply Delta approximation, 0= No delta
    inputs.solzen = (
        params.solzen
    )  # if DIRECT give solar zenith angle between 0 and 89 degrees (from 0 = nadir, 90 = horizon)

    # CHOOSE ATMOSPHERIC PROFILE for surface-incident flux:
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
    # For granular layers only, choose TOON
    # For granular layers + Fresnel layers below, choose ADD_DOUBLE
    ###############################################################

    inputs.TOON = False  # toggle Toon et al tridiagonal matrix solver
    inputs.ADD_DOUBLE = True  # toggle adding-doubling solver

    inputs.dz = params.dz  # thickness of each vertical layer (unit = m)
    inputs.nbr_lyr = len(params.dz)  # number of snow layers
    inputs.layer_type = (
        params.layer_type
    )  # Fresnel layers for the ADD_DOUBLE option, set all to 0 for the TOON option
    inputs.rho_layers = params.rho_layers  # density of each layer (unit = kg m-3)
    inputs.nbr_wvl = 480
    inputs.R_sfc = np.array(
        [0.25 for i in range(inputs.nbr_wvl)]
    )  # reflectance of undrlying surface - set across all wavelengths
    # inputs.R_sfc = np.genfromtxt('/home/joe/Code/BioSNICAR_GO_PY/Data/rain_polished_ice_spectrum.csv', delimiter = 'csv')

    ###############################################################################
    ## 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
    # For hexagonal plates or columns of any size choose GeometricOptics
    # For sphere, spheroids, koch snowflake with optional water coating choose Mie
    ###############################################################################

    inputs.rf_ice = 2  # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016

    # Ice grain shape can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3 = koch snowflake, 4 = hexagonal prisms
    # For 0,1,2,3:
    inputs.grain_shp = [0] * len(params.dz)  # grain shape(He et al. 2016, 2017)
    inputs.grain_rds = params.grain_rds  # effective grain radius of snow/bubbly ice
    inputs.rwater = [0] * len(params.dz)  # radius of optional liquid water coating

    # For 4:
    inputs.side_length = 0
    inputs.depth = 0

    # Shape factor = ratio of nonspherical grain effective radii to that of equal-volume sphere
    ### only activated when sno_shp > 1 (i.e. nonspherical)
    ### 0=use recommended default value (He et al. 2017)
    ### use user-specified value (between 0 and 1)
    inputs.shp_fctr = [0] * len(params.dz)

    # Aspect ratio (ratio of width to length)
    inputs.grain_ar = [0] * len(params.dz)

    #######################################
    ## 5) SET LAP CHARACTERISTICS
    #######################################

    # Define total number of different LAPs/aerosols in model
    inputs.nbr_aer = 30

    inputs.GA_units = 0
    inputs.SA_units = 0

    # determine C_factor (can be None or a number)
    # this is the concentrating factor that accounts for
    # resolution difference in field samples and model layers
    inputs.Cfactor_SA = 10
    inputs.Cfactor_GA = 10

    inputs.cdom_layer = [0] * len(params.dz)

    # Set names of files containing the optical properties of these LAPs:
    # inputs.FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
    inputs.FILE_soot1 = "bc_ChCB_rn40_dns1270.nc"
    inputs.FILE_soot2 = "miecot_slfsot_ChC90_dns_1317.nc"
    inputs.FILE_brwnC1 = "brC_Kirch_BCsd.nc"
    inputs.FILE_brwnC2 = "brC_Kirch_BCsd_slfcot.nc"
    inputs.FILE_dust1 = "dust_balkanski_central_size1.nc"
    inputs.FILE_dust2 = "dust_balkanski_central_size2.nc"
    inputs.FILE_dust3 = "dust_balkanski_central_size3.nc"
    inputs.FILE_dust4 = "dust_balkanski_central_size4.nc"
    inputs.FILE_dust5 = "dust_balkanski_central_size5.nc"
    inputs.FILE_ash1 = "volc_ash_eyja_central_size1.nc"
    inputs.FILE_ash2 = "volc_ash_eyja_central_size2.nc"
    inputs.FILE_ash3 = "volc_ash_eyja_central_size3.nc"
    inputs.FILE_ash4 = "volc_ash_eyja_central_size4.nc"
    inputs.FILE_ash5 = "volc_ash_eyja_central_size5.nc"
    inputs.FILE_ash_st_helens = "volc_ash_mtsthelens_20081011.nc"
    inputs.FILE_Skiles_dust1 = "dust_skiles_size1.nc"
    inputs.FILE_Skiles_dust2 = "dust_skiles_size2.nc"
    inputs.FILE_Skiles_dust3 = "dust_skiles_size3.nc"
    inputs.FILE_Skiles_dust4 = "dust_skiles_size4.nc"
    inputs.FILE_Skiles_dust5 = "dust_skiles_size5.nc"
    inputs.FILE_GreenlandCentral1 = "dust_greenland_central_size1.nc"
    inputs.FILE_GreenlandCentral2 = "dust_greenland_central_size2.nc"
    inputs.FILE_GreenlandCentral3 = "dust_greenland_central_size3.nc"
    inputs.FILE_GreenlandCentral4 = "dust_greenland_central_size4.nc"
    inputs.FILE_GreenlandCentral5 = "dust_greenland_central_size5.nc"
    inputs.FILE_Cook_Greenland_dust_L = "dust_greenland_Cook_LOW_20190911.nc"
    inputs.FILE_Cook_Greenland_dust_C = "dust_greenland_Cook_CENTRAL_20190911.nc"
    inputs.FILE_Cook_Greenland_dust_H = "dust_greenland_Cook_HIGH_20190911.nc"
    inputs.FILE_snw_alg = "snw_alg_r025um_chla020_chlb025_cara150_carb140.nc"
    inputs.FILE_glacier_algae = "Cook2020_glacier_algae_4_40.nc"

    # Add more glacier algae (not functional in current code)
    # (optical properties generated with GO), not included in the current model
    # algae1_r = 6 # algae radius
    # algae1_l = 60 # algae length
    # FILE_glacier_algae1 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae1_r,algae1_l))
    # algae2_r = 2 # algae radius
    # algae2_l = 10 # algae length
    # FILE_glacier_algae2 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae2_r,algae2_l))

    inputs.mss_cnc_soot1 = (
        params.mss_cnc_soot1
    )  # uncoated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_soot2 = [0] * len(
        params.dz
    )  # coated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_brwnC1 = [0] * len(
        params.dz
    )  # uncoated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_brwnC2 = [0] * len(
        params.dz
    )  # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_dust1 = (
        params.mss_cnc_dust1
    )  # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    inputs.mss_cnc_dust2 = [0] * len(
        params.dz
    )  # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    inputs.mss_cnc_dust3 = [0] * len(
        params.dz
    )  # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    inputs.mss_cnc_dust4 = [0] * len(
        params.dz
    )  # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    inputs.mss_cnc_dust5 = (
        params.mss_cnc_dust1
    )  # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    inputs.mss_cnc_ash1 = [0] * len(
        params.dz
    )  # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    inputs.mss_cnc_ash2 = [0] * len(
        params.dz
    )  # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    inputs.mss_cnc_ash3 = [0] * len(
        params.dz
    )  # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    inputs.mss_cnc_ash4 = [0] * len(
        params.dz
    )  # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    inputs.mss_cnc_ash5 = [0] * len(
        params.dz
    )  # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    inputs.mss_cnc_ash_st_helens = [0] * len(params.dz)  # ash from Mount Saint Helen's
    inputs.mss_cnc_Skiles_dust1 = [0] * len(
        params.dz
    )  # Colorado dust size 1 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust2 = [0] * len(
        params.dz
    )  # Colorado dust size 2 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust3 = [0] * len(
        params.dz
    )  # Colorado dust size 3 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust4 = [0] * len(
        params.dz
    )  # Colorado dust size 4 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust5 = [0] * len(
        params.dz
    )  # Colorado dust size 5 (Skiles et al 2017)
    inputs.mss_cnc_GreenlandCentral1 = [0] * len(
        params.dz
    )  # Greenland Central dust size 1 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral2 = [0] * len(
        params.dz
    )  # Greenland Central dust size 2 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral3 = [0] * len(
        params.dz
    )  # Greenland Central dust size 3 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral4 = [0] * len(
        params.dz
    )  # Greenland Central dust size 4 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral5 = [0] * len(
        params.dz
    )  # Greenland Central dust size 5 (Polashenski et al 2015)
    inputs.mss_cnc_Cook_Greenland_dust_L = [0] * len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_C = [0] * len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_H = [0] * len(params.dz)
    inputs.mss_cnc_snw_alg = [0] * len(
        params.dz
    )  # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
    inputs.mss_cnc_glacier_algae = (
        params.mss_cnc_glacier_algae
    )  # glacier algae type1 (Cook et al. 2020)

    outputs = snicar_feeder(inputs)

    return outputs.albedo, outputs.BBA
