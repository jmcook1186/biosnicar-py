import collections as c
import sys

sys.path.append("./src")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from snicar_feeder import snicar_feeder


def generate_snicar_params(
    layer_type, density, dz, alg, solzen, reff, bc, dust1, dust5
):

    rho = [density] * len(dz)
    reff = [reff] * len(dz)
    layer_type = [layer_type] * len(dz)
    mss_cnc_glacier_algae = [0] * len(dz)
    mss_cnc_soot1 = [bc]*len(dz)
    mss_cnc_dust1 = [0] * len(dz)
    mss_cnc_dust1[0] = dust1
    mss_cnc_dust5 = [0] * len(dz)
    mss_cnc_dust5[0] = dust5

    params = c.namedtuple(
        "params",
        "rho, shp, lyr, toon, add_double, incoming, grain_rds, direct, aprx_typ, layer_type, dz, mss_cnc_glacier_algae, bc, mss_cnc_dust1, mss_cnc_dust5, solzen",
    )
    params.lyr = 0
    params.rds = reff
    params.rho = rho
    params.layer_type = layer_type
    params.dz = dz
    params.mss_cnc_glacier_algae = mss_cnc_glacier_algae
    params.solzen = solzen
    params.bc = mss_cnc_soot1
    params.mss_cnc_dust1 = mss_cnc_dust1
    params.mss_cnc_dust5 = mss_cnc_dust5
    params.direct = 1
    params.aprx_typ = 3
    params.incoming = 4
    params.toon = 0
    params.add_double = 1
    params.shp = 0
    return params


def call_snicar(params):
    
    # --------------------------------------------------------------------------------------
    # 1) Initialize Inputs of the model
    # --------------------------------------------------------------------------------------

    inputs = c.namedtuple(
        "Inputs",
        [
            "dir_base",
            "verbosity",
            "rf_ice",
            "mie_solver",
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
            "lwc_pct_bubbles",
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
    dir_base_path = Path(__file__).parent.parent
    inputs.dir_base = str(dir_base_path) + "/"
    savepath = inputs.dir_base  # base path for saving figures
    WRITE_CONFIG_TO_TEXTFILE = False  # toggle to save config to file
    inputs.verbosity = 0  # 1 to print real-time updates


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

    inputs.direct = (
        params.direct
    )  # 1= direct-beam incident flux, 0= Diffuse incident flux
    inputs.aprx_typ = (
        params.aprx_typ
    )  # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    inputs.delta = 1  # 1= Apply delta approximation, 0= No delta
    inputs.solzen = (
        params.solzen
    )  # if direct give solar zenith angle between 0 and 89 degrees (from 0 = nadir, 90 = horizon)

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
    inputs.incoming_i = params.incoming
    
    # --------------------------------------------------------------------------------------
    # 4) SET UP ICE/SNOW LAYERS
    # For granular layers only, choose toon
    # For granular layers + Fresnel layers, choose add_double
    # --------------------------------------------------------------------------------------

    inputs.toon = params.toon  # toggle toon et al tridiagonal matrix solver
    inputs.add_double = params.add_double  # toggle adding-doubling solver

    inputs.dz = params.dz  # thickness of each vertical layer (unit = m)
    inputs.nbr_lyr = len(params.dz)  # number of snow layers
    inputs.layer_type = [params.lyr] * len(
        params.dz
    )  
    # 0 = ice grain layer of different shapes,
    # 1 = solid bubbly ice w/ fresnel layer on top,
    # 2 = liquid water film,
    # 3 = solid bubbly ice w/out fresnel layer
    # 4 = wet snow mixing ice and water spheres
    inputs.lwc = [0] * len(inputs.dz)
    inputs.lwc_pct_bubbles = 0  # amount of water as bubbles

    inputs.cdom_layer = [0] * len(inputs.dz)  # Only for layer type == 1
    inputs.rho_layers = params.rho  # density of each layer (unit = kg m-3)
    inputs.nbr_wvl = 480

    # reflectance of underlying surface - set across all wavelengths
    rain_polished_ice_spectrum_file = (
        dir_base_path.joinpath("Data")
        .joinpath("OP_data")
        .joinpath("480band")
        .joinpath("r_sfc")
        .joinpath("rain_polished_ice_spectrum.csv")
    )

    inputs.R_sfc = np.genfromtxt(
        str(rain_polished_ice_spectrum_file),
        delimiter="csv",
    )  # import underlying ice from file


    # --------------------------------------------------------------------------------------
    # 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
    # For hexagonal plates or columns of any size choose GeometricOptics
    # For sphere, spheroids, koch snowflake with optional water coating choose Mie
    # --------------------------------------------------------------------------------------

    # define source of ice refractive index data.
    # 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
    inputs.rf_ice = 2
    
    # define solver for Mie calculations
    # 0 = Wiscombe 1979 (default)
    # 1 = Bohren and Huffman 1983 (as per the Matlab version)
    inputs.mie_solver = 1

    # Ice grain shape can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3 = koch snowflake, 4 = hexagonal prisms
    # For 0,1,2,3:
    inputs.grain_shp = [params.shp] * len(
        params.dz
    )  # grain shape(He et al. 2016, 2017)
    inputs.grain_rds = params.rds  # effective grain radius of snow/bubbly ice
    inputs.rwater = [0] * len(params.dz)  # radius of optional liquid water coating

    # For 4:
    inputs.side_length = [10000, 10000]
    inputs.depth = [10000, 10000]

    # Shape factor = ratio of nonspherical grain effective radii to that of equal-volume sphere
    ### only activated when sno_shp > 1 (i.e. nonspherical)
    ### 0=use recommended default value (He et al. 2017)
    ### use user-specified value (between 0 and 1)
    inputs.shp_fctr = [0] * len(params.dz)

    # Aspect ratio (ratio of width to length)
    inputs.grain_ar = [0] * len(params.dz)

    # --------------------------------------------------------------------------------------
    # 5) SET LAP CHARACTERISTICS
    # --------------------------------------------------------------------------------------

    # Define total number of different LAPs/aerosols in model
    inputs.nbr_aer = 30

    # define units for algae absorption cross section input file
    # 0 = m2/kg for MAC, ppb for mss_cnc (this is default)
    # 1 = m2/cell for MAC, cell/mL for mss_cnc
    inputs.GA_units = 1  # glacier algae
    inputs.SA_units = 1  # snow algae

    # Set names of files containing the optical properties of these LAPs:
    inputs.file_soot1 = (
        "bc_ChCB_rn40_dns1270.nc"  # uncoated BC (Bohren and Huffman, 1983) 
    )
    inputs.file_soot2 = (
        "bc_ChCB_rn40_dns1270_slfcot.nc"  # coated BC (Bohren and Huffman, 1983)
    ) 
    inputs.file_brwnC1 = (
        "brC_Kirch_BCsd.nc"  # uncoated brown carbon (Kirchstetter et al. (2004).)
    )
    inputs.file_brwnC2 = "brC_Kirch_BCsd_slfcot.nc"  # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.file_dust1 = "dust_balkanski_central_size1.nc"  # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    inputs.file_dust2 = "dust_balkanski_central_size2.nc"  # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    inputs.file_dust3 = "dust_balkanski_central_size3.nc"  # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    inputs.file_dust4 = "dust_balkanski_central_size4.nc"  # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    inputs.file_dust5 = "dust_balkanski_central_size5.nc"  # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    inputs.file_ash1 = "volc_ash_eyja_central_size1.nc"  # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    inputs.file_ash2 = "volc_ash_eyja_central_size2.nc"  # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    inputs.file_ash3 = "volc_ash_eyja_central_size3.nc"  # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    inputs.file_ash4 = "volc_ash_eyja_central_size4.nc"  # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    inputs.file_ash5 = "volc_ash_eyja_central_size5.nc"  # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    inputs.file_ash_st_helens = (
        "volc_ash_mtsthelens_20081011.nc"  # ashes from Mount Saint Helen's
    )
    inputs.file_Skiles_dust1 = (
        "dust_skiles_size1.nc"  # Colorado dust 1 (Skiles et al 2017)
    )
    inputs.file_Skiles_dust2 = (
        "dust_skiles_size2.nc"  # Colorado dust 2 (Skiles et al 2017)
    )
    inputs.file_Skiles_dust3 = (
        "dust_skiles_size3.nc"  # Colorado dust 3 (Skiles et al 2017)
    )
    inputs.file_Skiles_dust4 = (
        "dust_skiles_size4.nc"  # Colorado dust 4 (Skiles et al 2017)
    )
    inputs.file_Skiles_dust5 = (
        "dust_skiles_size5.nc"  # Colorado dust 5 (Skiles et al 2017)
    )
    inputs.file_GreenlandCentral1 = (
        "dust_greenland_central_size1.nc"  # Greenland dust 1 (Polashenski et al 2015)
    )
    inputs.file_GreenlandCentral2 = (
        "dust_greenland_central_size2.nc"  # Greenland dust 2 (Polashenski et al 2015)
    )
    inputs.file_GreenlandCentral3 = (
        "dust_greenland_central_size3.nc"  # Greenland dust 3 (Polashenski et al 2015)
    )
    inputs.file_GreenlandCentral4 = (
        "dust_greenland_central_size4.nc"  # Greenland dust 4 (Polashenski et al 2015)
    )
    inputs.file_GreenlandCentral5 = (
        "dust_greenland_central_size5.nc"  # Greenlanddust 5 (Polashenski et al 2015)
    )
    inputs.file_Cook_Greenland_dust_L = "dust_greenland_Cook_LOW_20190911.nc"  # GRIS dust (Cook et al. 2019 "LOW") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.file_Cook_Greenland_dust_C = "dust_greenland_Cook_CENTRAL_20190911.nc"  # GRIS dust 1 (Cook et al. 2019 "mean") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.file_Cook_Greenland_dust_H = "dust_greenland_Cook_HIGH_20190911.nc"  # GRIS dust 1 (Cook et al. 2019 "HIGH") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.file_snw_alg = "snow_algae_empirical_Chevrollier2023.nc"  # Snow Algae (spherical, C nivalis)
    inputs.file_glacier_algae = "ice_algae_empirical_Chevrollier2023.nc"  # glacier algae in cells/ml or ppb depending on GA_units

    # Add more glacier algae (not functional in current code)
    # (optical properties generated with GO), not included in the current model
    # algae1_r = 6 # algae radius
    # algae1_l = 60 # algae length
    # file_glacier_algae1 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae1_r,algae1_l))
    # algae2_r = 2 # algae radius
    # algae2_l = 10 # algae length
    # file_glacier_algae2 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae2_r,algae2_l))

    inputs.mss_cnc_soot1 = params.bc  # uncoated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_soot2 = [0] * len(
        params.dz
    )  # coated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_brwnC1 = [0] * len(
        params.dz
    )  # uncoated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_brwnC2 = [0] * len(
        params.dz
    )  # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_dust1 = [0] * len(
        params.dz
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
    inputs.mss_cnc_dust5 = [0] * len(
        params.dz
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
    inputs.mss_cnc_glacier_algae = [0] * len(
        params.dz
    )  # glacier algae type1 (Cook et al. 2020)

    outputs = snicar_feeder(inputs)

    return outputs.albedo, outputs.BBA
