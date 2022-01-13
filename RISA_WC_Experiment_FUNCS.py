"""
Joseph Cook, Aarhus University, Feb 2021

This script contains functions that a) generate 
SNICAR-predicted spectral albedo that approximate
field-measured spectra for a variety of weathering 
crust configurations; b) quantify the albedo change
resulting from a range of WC development scenarios 


includes:

find_best_params()
    the forward modelling script that retrieves the snicar params that generate the
    best-matching curve to a given field spectrum

call_snicar()
    the function used to do multiple snicar runs with params provided as a named tuple

match_field_spectra()
    the function for plotting simulated and measured field spectra in a multipanel fig

isolate_biological_effect()
    the function for calculating the albedo reduction due to bio vs phys processes by spectral differencing

build_LUT()
    function for constructing lookup table to be used in the inverse model

inverse_model()
    function finds the best matching entry in the LUTs for each field spectrum and returns the params
    used to generate it


"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SNICAR_feeder import snicar_feeder
import statsmodels.api as sm  
import collections as c
import xarray as xr
import dask


def find_best_params(field_data_fname, sampleID, CIsites, LAsites, HAsites, weight, clean=True):

    """
    This function will return the SNICAR parameter set that provides the
    best approximation to a given field spectrum. 

    """

    spectra = pd.read_csv(field_data_fname)
    spectra = spectra[::10]

    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)]

    if sampleID =='CImean':
        field_spectrum = CIspec.mean(axis=1)
    elif sampleID =='HAmean':
        field_spectrum = HAspec.mean(axis=1)
    elif sampleID == 'LAmean':
        field_spectrum = LAspec.mean(axis=1)
    elif sampleID == 'RAIN':
        field_spectrum = spectra['RAIN2']
    else:
        field_spectrum = spectra[sampleID]

    # calculate 2BDA index of field spectrum
    BDA2idx = np.array(field_spectrum)[36]/np.array(field_spectrum)[31]

    dens = [550, 600, 650, 700, 750, 800, 850]
    dz = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2]
    rds = [600, 700, 800, 900, 1000]
    alg = [0, 2500, 5000, 7500, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000]
    Cfactors = [10, 20, 30]
    solzen = [40, 45, 50]

    @dask.delayed
    def run_sims(i,j,k,p,q,z):
        
        params = c.namedtuple("params","rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, Cfactor_GA, solzen")
        params.rho_layers = [i,i]
        params.grain_rds = [j,j]
        params.layer_type = [1,1]
        params.dz = [0.001,k]
        params.mss_cnc_glacier_algae = [p,0]
        params.Cfactor_GA = Cfactors[q]
        params.solzen = z

        albedo, BBA = call_snicar(params)
        
        error_vis = np.mean(abs(albedo[15:55]-field_spectrum[0:40]))
        error_nir = np.mean(abs(albedo[55:100]-field_spectrum[40:85]))
        error = ((error_vis*weight)+error_nir)/(1+weight)

        params = (i,j,k,p,z)
        out =(params, error)

        return out

    Out = []
    # now use the reduced LUT to call snicar and obtain best matching spectrum
    for i in dens:
        for j in rds:
            for k in dz:
                for p in alg:
                    for z in solzen:
                        out = run_sims(i,j,k,p,z)
                        Out.append(out)
    
    Result = dask.compute(*Out,scheduler='processes')

    return Result


def call_snicar(params):
    
    inputs = c.namedtuple('inputs',['dir_base', 'verbosity',\
        'rf_ice', 'incoming_i', 'DIRECT', 'layer_type',\
        'cdom_layer', 'APRX_TYP', 'DELTA', 'solzen', \
        'TOON', 'ADD_DOUBLE', 'R_sfc', 'dz', 'rho_layers',\
        'grain_rds','side_length', 'depth', 'rwater', 'nbr_lyr',\
        'nbr_aer', 'grain_shp', 'shp_fctr', 'grain_ar', 'GA_units',\
        'SA_units', 'Cfactor_GA','Cfactor_SA','mss_cnc_soot1',\
        'mss_cnc_soot2', 'mss_cnc_brwnC1', 'mss_cnc_brwnC2',\
        'mss_cnc_dust1', 'mss_cnc_dust2', 'mss_cnc_dust3',\
        'mss_cnc_dust4', 'mss_cnc_dust5', 'mss_cnc_ash1',\
        'mss_cnc_ash2', 'mss_cnc_ash3', 'mss_cnc_ash4', \
        'mss_cnc_ash5', 'mss_cnc_ash_st_helens',\
        'mss_cnc_Skiles_dust1', 'mss_cnc_Skiles_dust2',\
        'mss_cnc_Skiles_dust3', 'mss_cnc_Skiles_dust4',\
        'mss_cnc_Skiles_dust5', 'mss_cnc_GreenlandCentral1',\
        'mss_cnc_GreenlandCentral2', 'mss_cnc_GreenlandCentral3',\
        'mss_cnc_GreenlandCentral4', 'mss_cnc_GreenlandCentral5',\
        'mss_cnc_Cook_Greenland_dust_L', 'mss_cnc_Cook_Greenland_dust_C',\
        'mss_cnc_Cook_Greenland_dust_H', 'mss_cnc_snw_alg',\
        'mss_cnc_glacier_algae', 'FILE_soot1', 'FILE_soot2',\
        'FILE_brwnC1', 'FILE_brwnC2', 'FILE_dust1', 'FILE_dust2',\
        'FILE_dust3', 'FILE_dust4', 'FILE_dust5', 'FILE_ash1',\
        'FILE_ash2', 'FILE_ash3', 'FILE_ash4', 'FILE_ash5',\
        'FILE_ash_st_helens', 'FILE_Skiles_dust1', 'FILE_Skiles_dust2',\
        'FILE_Skiles_dust3', 'FILE_Skiles_dust4', 'FILE_Skiles_dust5',\
        'FILE_GreenlandCentral1', 'FILE_GreenlandCentral2',\
        'FILE_GreenlandCentral3', 'FILE_GreenlandCentral4',\
        'FILE_GreenlandCentral5', 'FILE_Cook_Greenland_dust_L',\
        'FILE_Cook_Greenland_dust_C', 'FILE_Cook_Greenland_dust_H',\
        'FILE_snw_alg', 'FILE_glacier_algae', 'tau', 'g', 'SSA',\
        'mu_not', 'nbr_wvl', 'wvl', 'Fs', 'Fd', 'L_snw', 'flx_slr'])


    ##############################
    ## 2) Set working directory 
    ##############################

    # set dir_base to the location of the BioSNICAR_GO_PY folder
    inputs.dir_base = '/home/joe/Code/BioSNICAR_GO_PY/'
    savepath = inputs.dir_base # base path for saving figures
    write_config_to_textfile = False # toggle to save config to file
    inputs.verbosity = 0 # 1 to print real-time updates

    ################################
    ## 3) Choose plot/print options
    ################################

    show_figs = True # toggle to display spectral albedo figure
    save_figs = False # toggle to save spectral albedo figure to file
    print_BBA = True # toggle to print broadband albedo to terminal
    print_band_ratios = False # toggle to print various band ratios to terminal
    smooth = False # apply optional smoothing function (Savitzky-Golay filter)
    window_size = 9 # if applying smoothing filter, define window size
    poly_order = 3 # if applying smoothing filter, define order of polynomial

    #######################################
    ## 4) RADIATIVE TRANSFER CONFIGURATION
    #######################################

    inputs.DIRECT   = 1       # 1= Direct-beam incident flux, 0= Diffuse incident flux
    inputs.APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    inputs.DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
    inputs.solzen   = params.solzen      # if DIRECT give solar zenith angle between 0 and 89 degrees (from 0 = nadir, 90 = horizon)

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

    inputs.TOON = False # toggle Toon et al tridiagonal matrix solver
    inputs.ADD_DOUBLE = True # toggle adding-doubling solver

    inputs.dz = params.dz # thickness of each vertical layer (unit = m)
    inputs.nbr_lyr = len(params.dz)  # number of snow layers
    inputs.layer_type = params.layer_type # Fresnel layers for the ADD_DOUBLE option, set all to 0 for the TOON option
    inputs.cdom_layer = [0,0] # Only for layer type == 1, CDOM data from L Halbach
    inputs.rho_layers = params.rho_layers # density of each layer (unit = kg m-3) 
    inputs.nbr_wvl=480 
    #inputs.R_sfc = np.array([0.1 for i in range(inputs.nbr_wvl)]) # reflectance of undrlying surface - set across all wavelengths
    inputs.R_sfc = np.genfromtxt('./Data/rain_polished_ice_spectrum.csv', delimiter = 'csv') # import underlying ice from file


    ###############################################################################
    ## 5) SET UP OPTICAL & PHYSICAL PROPERTIES OF SNOW/ICE GRAINS
    # For hexagonal plates or columns of any size choose GeometricOptics
    # For sphere, spheroids, koch snowflake with optional water coating choose Mie
    ###############################################################################

    inputs.rf_ice = 2 # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016

    # Ice grain shape can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3 = koch snowflake, 4 = hexagonal prisms
    # For 0,1,2,3:
    inputs.grain_shp =[0]*len(params.dz) # grain shape(He et al. 2016, 2017)
    inputs.grain_rds = params.grain_rds # effective grain radius of snow/bubbly ice
    inputs.rwater = [0]*len(params.dz) # radius of optional liquid water coating
    
    # For 4:
    inputs.side_length = 0 
    inputs.depth = 0

    # Shape factor = ratio of nonspherical grain effective radii to that of equal-volume sphere
    ### only activated when sno_shp > 1 (i.e. nonspherical)
    ### 0=use recommended default value (He et al. 2017)
    ### use user-specified value (between 0 and 1)
    inputs.shp_fctr = [0]*len(params.dz) 

    # Aspect ratio (ratio of width to length)
    inputs.grain_ar = [0]*len(params.dz) 

    #######################################
    ## 5) SET LAP CHARACTERISTICS
    #######################################

    # Define total number of different LAPs/aerosols in model
    inputs.nbr_aer = 30

    # define units for algae absorption cross section input file
    # 0 = m2/kg for MAC, ppb for mss_cnc (this is default)
    # 1 = m2/cell for MAC, cell/mL for mss_cnc
    inputs.GA_units = 1 # glacier algae
    inputs.SA_units = 1 # snow algae

    # determine C_factor (can be None or a number)
    # this is the concentrating factor that accounts for
    # resolution difference in field samples and model layers
    inputs.Cfactor_GA = params.Cfactor_GA
    inputs.Cfactor_SA = 0

    # Set names of files containing the optical properties of these LAPs:
    inputs.FILE_soot1 ='mie_sot_ChC90_dns_1317.nc' # uncoated BC (Bohren and Huffman, 1983)
    inputs.FILE_soot2 = 'miecot_slfsot_ChC90_dns_1317.nc' # coated BC (Bohren and Huffman, 1983)
    inputs.FILE_brwnC1 = 'brC_Kirch_BCsd.nc' # uncoated brown carbon (Kirchstetter et al. (2004).)
    inputs.FILE_brwnC2 = 'brC_Kirch_BCsd_slfcot.nc' # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.FILE_dust1 = 'dust_balkanski_central_size1.nc' # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    inputs.FILE_dust2 = 'dust_balkanski_central_size2.nc' # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    inputs.FILE_dust3 = 'dust_balkanski_central_size3.nc' # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    inputs.FILE_dust4 = 'dust_balkanski_central_size4.nc' # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    inputs.FILE_dust5 = 'dust_balkanski_central_size5.nc' # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    inputs.FILE_ash1 = 'volc_ash_eyja_central_size1.nc' # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    inputs.FILE_ash2 = 'volc_ash_eyja_central_size2.nc' # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    inputs.FILE_ash3 = 'volc_ash_eyja_central_size3.nc' # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    inputs.FILE_ash4 = 'volc_ash_eyja_central_size4.nc' # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    inputs.FILE_ash5 = 'volc_ash_eyja_central_size5.nc' # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    inputs.FILE_ash_st_helens = 'volc_ash_mtsthelens_20081011.nc' # ashes from Mount Saint Helen's
    inputs.FILE_Skiles_dust1 = 'dust_skiles_size1.nc' # Colorado dust 1 (Skiles et al 2017)
    inputs.FILE_Skiles_dust2 = 'dust_skiles_size2.nc' # Colorado dust 2 (Skiles et al 2017)
    inputs.FILE_Skiles_dust3 = 'dust_skiles_size3.nc' # Colorado dust 3 (Skiles et al 2017)
    inputs.FILE_Skiles_dust4 = 'dust_skiles_size4.nc' # Colorado dust 4 (Skiles et al 2017)
    inputs.FILE_Skiles_dust5 = 'dust_skiles_size5.nc' # Colorado dust 5 (Skiles et al 2017)
    inputs.FILE_GreenlandCentral1 = 'dust_greenland_central_size1.nc' # Greenland dust 1 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral2 = 'dust_greenland_central_size2.nc' # Greenland dust 2 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral3 = 'dust_greenland_central_size3.nc' # Greenland dust 3 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral4 = 'dust_greenland_central_size4.nc' # Greenland dust 4 (Polashenski et al 2015)
    inputs.FILE_GreenlandCentral5  = 'dust_greenland_central_size5.nc' # Greenlanddust 5 (Polashenski et al 2015)
    inputs.FILE_Cook_Greenland_dust_L = 'dust_greenland_Cook_LOW_20190911.nc' # GRIS dust (Cook et al. 2019 "LOW") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.FILE_Cook_Greenland_dust_C = 'dust_greenland_Cook_CENTRAL_20190911.nc' # GRIS dust 1 (Cook et al. 2019 "mean") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.FILE_Cook_Greenland_dust_H = 'dust_greenland_Cook_HIGH_20190911.nc' # GRIS dust 1 (Cook et al. 2019 "HIGH") NOT FUNCTIONAL IN THIS RELEASE (COMING SOON)
    inputs.FILE_snw_alg  = 'SA_Halbach2021_oct21.nc' # Snow Algae (spherical, C nivalis)
    inputs.FILE_glacier_algae = 'GA_Chevrollier2022_r4.9_L18.8.nc' # glacier algae in cells/ml or ppb depending on GA_units 

    # Add more glacier algae (not functional in current code)
    # (optical properties generated with GO), not included in the current model
    # algae1_r = 6 # algae radius
    # algae1_l = 60 # algae length
    # FILE_glacier_algae1 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae1_r,algae1_l))
    # algae2_r = 2 # algae radius
    # algae2_l = 10 # algae length
    # FILE_glacier_algae2 = str(dir_go_lap_files + 'RealPhenol_algae_geom_{}_{}.nc'.format(algae2_r,algae2_l))

        
    inputs.mss_cnc_soot1 = [0]*len(params.dz)    # uncoated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_soot2 = [0]*len(params.dz)    # coated black carbon (Bohren and Huffman, 1983)
    inputs.mss_cnc_brwnC1 = [0]*len(params.dz)   # uncoated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_brwnC2 = [0]*len(params.dz)   # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    inputs.mss_cnc_dust1 = [0]*len(params.dz)    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    inputs.mss_cnc_dust2 = [0]*len(params.dz)    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    inputs.mss_cnc_dust3 = [0]*len(params.dz)    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    inputs.mss_cnc_dust4 = [0]*len(params.dz)    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    inputs.mss_cnc_dust5 = [0]*len(params.dz)    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    inputs.mss_cnc_ash1 = [0]*len(params.dz)    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    inputs.mss_cnc_ash2 = [0]*len(params.dz)    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    inputs.mss_cnc_ash3 = [0]*len(params.dz)    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    inputs.mss_cnc_ash4 = [0]*len(params.dz)    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    inputs.mss_cnc_ash5 = [0]*len(params.dz)    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    inputs.mss_cnc_ash_st_helens = [0]*len(params.dz)   # ash from Mount Saint Helen's
    inputs.mss_cnc_Skiles_dust1 = [0]*len(params.dz)    # Colorado dust size 1 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust2 = [0]*len(params.dz)    # Colorado dust size 2 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust3 = [0]*len(params.dz)    # Colorado dust size 3 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust4 = [0]*len(params.dz)  # Colorado dust size 4 (Skiles et al 2017)
    inputs.mss_cnc_Skiles_dust5 = [0]*len(params.dz)  # Colorado dust size 5 (Skiles et al 2017)
    inputs.mss_cnc_GreenlandCentral1 = [0]*len(params.dz) # Greenland Central dust size 1 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral2 = [0]*len(params.dz) # Greenland Central dust size 2 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral3 = [0]*len(params.dz) # Greenland Central dust size 3 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral4 = [0]*len(params.dz) # Greenland Central dust size 4 (Polashenski et al 2015)
    inputs.mss_cnc_GreenlandCentral5 = [0]*len(params.dz) # Greenland Central dust size 5 (Polashenski et al 2015)
    inputs.mss_cnc_Cook_Greenland_dust_L = [0]*len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_C = [0]*len(params.dz)
    inputs.mss_cnc_Cook_Greenland_dust_H = [0]*len(params.dz)
    inputs.mss_cnc_snw_alg = [0]*len(params.dz)    # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
    inputs.mss_cnc_glacier_algae = params.mss_cnc_glacier_algae   # glacier algae type1 (Cook et al. 2020)

    nbr_aer = 30
    
    outputs = snicar_feeder(inputs)
    

    return outputs.albedo, outputs.BBA


def match_field_spectra(field_data_fname, fnames, rho, rds, dz, alg, measured_cells,\
    CIsites, LAsites, HAsites, apply_ARF, plot_ARF, ARF_CI, ARF_HA, savepath):
    
    """
    plot field against SNICAR spectra
    requires parameters to be known in advance and hard coded inside this function
    the relevant params can be generated using the find_best_params() func

    params:
    field_data_fname = filename for spectral database
    
    The following params are parallel arrays - the order matters and 
    must match the filenames! Pairs of values represent values for 
    upper and lower layers in model. These are the values used to
    generate snicar spectrum to match field spectrum with name = fname[i]
    
    fnames = sample IDs for spectra to match model runs
    rho = pairs of density values [upper, lower]
    rds = pairs of r_eff values [upper, lower]
    dz = layer thickness [upper, lower] NB. upper is always 0.001
    alg = mass concentration of algae [upper, lower] NB. lower is always 0
    measured_cells = array of the actual measured cell concentration for each spectrum

    e.g.
    fnames= ['2016_WI_8','14_7_SB6','14_7_SB9','14_7_SB1','21_7_SB2','14_7_SB2', '22_7_SB3', 'RAIN']
    rho = [[550,550],[650,650],[800,800],[850,850],[750,750],[800,800],[800,800],[900,900]]
    rds = [[550,550],[650,650],[850,850],[850,850],[800,800],[800,800],[750,750],[900,900]]
    dz = [[0.001,0.3],[0.001,0.09],[0.001,0.03],[0.001,0.03],[0.001,0.02],[0.001,0.06],[0.001,0.05],[0.001,0.03]]
    alg = [[0,0],[0,0],[20000,0],[30000,0],[45000,0],[3000,0],[8000,0],[0,0]]  

    returns:
    None, but saves figure to savepath

    """

    spectra = pd.read_csv(field_data_fname)

    # reformat feld spectra to match snicar resolution
    spectra = spectra[::10]

    # gather spectra for each surface type
    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)]
    RAINspec = spectra['RAIN2']


    if plot_ARF:
        plt.plot(spectra.Wavelength[0:100],ARF_CI[0:100],marker='x', label='clean ice ARF'),
        plt.plot(spectra.Wavelength[0:100],ARF_HA[0:100],marker='o',linestyle='dashed',label='algal ice ARF')
        plt.ylabel('Anitostropic Reflectance Factor'),plt.xlabel("Wavelength (nm)")
        plt.legend(loc='best')
        plt.savefig(str(savepath+'ARF.jpg'),dpi=300)

    # define local function for calling snicar
    def simulate_albedo(rds, rho, dz, alg):
  
        params = c.namedtuple("params",\
            "rho_layers, grain_rds, layer_type, dz,\
                 mss_cnc_glacier_algae, Cfactor_GA, solzen")
        params.grain_rds = rds
        params.rho_layers = rho
        params.layer_type = [1,1]
        params.dz = dz
        params.mss_cnc_glacier_algae = alg
        params.Cfactor_GA = 30
        params.solzen = 40
        albedo, BBA = call_snicar(params)
  
        return albedo, BBA
  
   
    # set up output array
    # and call snicar with each set of params
    OutArray = np.zeros(shape=(len(fnames),480))

    for i in np.arange(0,len(fnames),1):

        albedo,BBA=simulate_albedo(rds[i],rho[i],dz[i],alg[i])
        if apply_ARF:
            if alg[i][0] > 5000:
                albedo[15:230] = albedo[15:230]*ARF_HA
            else:
                albedo[15:230] = albedo[15:230]*ARF_CI

        OutArray[i,:] = albedo


    # calculate mean absolute error for model vs measured spectrum
    error = []
    for i in np.arange(0,len(fnames),1):
        error.append(abs(np.mean(spectra[fnames[i]].iloc[0:130] - (OutArray[i,15:145]))))

    print("wavelength = ",spectra.Wavelength.iloc[130])

    # plot figure
    fig,axes = plt.subplots(4,2,figsize=(10,10))

    if apply_ARF:
        ylabel='Reflectance'
        PlotName = 'FieldvsMeasuredReflectance.jpg' 
    else:
        ylabel='Albedo'
        PlotName = 'FieldvsMeasuredAlbedo.jpg'

    axes[0,0].plot(spectra.Wavelength,spectra['23_7_SB1'],label='field')
    axes[0,0].plot(spectra.Wavelength,OutArray[0,15:230],label='model',linestyle='--')
    axes[0,0].set_ylim(0,1), axes[0,0].set_xlim(350,1800)
    axes[0,0].text(1400,0.2,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[0][0],rho[0][0],dz[0][1],error[0]))
    axes[0,0].text(400,0.1,'{} \n{} cells/mL '.format(fnames[0],measured_cells[0]))
    axes[0,0].set_ylabel(ylabel), axes[0,0].set_xlabel('Wavelength (nm)')
    axes[0,0].legend(loc='best')

    axes[0,1].plot(spectra.Wavelength,spectra['14_7_SB6'])
    axes[0,1].plot(spectra.Wavelength,OutArray[1,15:230],linestyle='--')
    axes[0,1].set_ylim(0,1), axes[0,1].set_xlim(350,1800)
    axes[0,1].text(1450,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[1][0],rho[1][0],dz[1][1],error[1]))
    axes[0,1].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[1],measured_cells[1]))
    axes[0,1].set_ylabel(ylabel), axes[0,1].set_xlabel('Wavelength (nm)')

    axes[1,0].plot(spectra.Wavelength,spectra['14_7_SB9'])
    axes[1,0].plot(spectra.Wavelength,OutArray[2,15:230],linestyle='--')
    axes[1,0].set_ylim(0,1), axes[1,0].set_xlim(350,1800)
    axes[1,0].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[2][0],rho[2][0],dz[2][1],error[2]))
    axes[1,0].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[2],measured_cells[2]))
    axes[1,0].set_ylabel(ylabel), axes[1,0].set_xlabel('Wavelength (nm)')

    axes[1,1].plot(spectra.Wavelength,spectra['14_7_SB1'])
    axes[1,1].plot(spectra.Wavelength,OutArray[3,15:230],linestyle='--')
    axes[1,1].set_ylim(0,1), axes[1,1].set_xlim(350,1800)
    axes[1,1].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[3][0],rho[3][0],dz[3][1],error[3]))
    axes[1,1].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[3],measured_cells[3]))
    axes[1,1].set_ylabel(ylabel), axes[1,1].set_xlabel('Wavelength (nm)')

    axes[2,0].plot(spectra.Wavelength,spectra['22_7_SB5'])
    axes[2,0].plot(spectra.Wavelength,OutArray[4,15:230],linestyle='--')
    axes[2,0].set_ylim(0,1), axes[2,0].set_xlim(350,1800)
    axes[2,0].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[4][0],rho[4][0],dz[4][1], error[4]))
    axes[2,0].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[4],measured_cells[4]))
    axes[2,0].set_ylabel(ylabel), axes[2,0].set_xlabel('Wavelength (nm)')

    axes[2,1].plot(spectra.Wavelength,spectra['14_7_SB2'])
    axes[2,1].plot(spectra.Wavelength,OutArray[5,15:230],linestyle='--')
    axes[2,1].set_ylim(0,1), axes[2,1].set_xlim(350,1800)
    axes[2,1].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[5][0],rho[5][0],dz[5][1],error[5]))
    axes[2,1].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[5],measured_cells[5]))
    axes[2,1].set_ylabel(ylabel), axes[2,1].set_xlabel('Wavelength (nm)')

    axes[3,0].plot(spectra.Wavelength,spectra['22_7_SB3'])
    axes[3,0].plot(spectra.Wavelength,OutArray[6,15:230],linestyle='--')
    axes[3,0].set_ylim(0,1), axes[3,0].set_xlim(350,1800)
    axes[3,0].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[6][0],rho[6][0],dz[6][1],error[6]))
    axes[3,0].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[6],measured_cells[6]))
    axes[3,0].set_ylabel(ylabel), axes[3,0].set_xlabel('Wavelength (nm)')

    axes[3,1].plot(spectra.Wavelength,RAINspec)
    axes[3,1].plot(spectra.Wavelength,OutArray[7,15:230],linestyle='--')
    axes[3,1].set_ylim(0,1), axes[3,1].set_xlim(350,1800)
    axes[3,1].text(1400,0.3,"r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(\
        rds[7][0],rho[7][0],dz[7][1],error[7]))
    axes[3,1].text(400,0.8,'{}: \n{} cells/mL'.format(fnames[1],measured_cells[1]))
    axes[3,1].set_ylabel(ylabel), axes[3,1].set_xlabel('Wavelength (nm)')

    fig.tight_layout()
    plt.savefig(str(savepath+PlotName),dpi=300)

    return



def isolate_biological_effect(field_data_fname, CIsites, LAsites, HAsites, savepath):

    """
    This function estimates the albedo reduction resulting from the ice physical changes
    versus the biological growth. 

    Some nuance to the interpretation because the ce surface likely would not
    degrade to the same extent without the algal bloom.

    """

    #read in spectral database
    spectra = pd.read_csv(field_data_fname)

    # reformat feld spectra to match snicar resolution
    spectra = spectra[::10]
    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)] 
    meanCI = CIspec.mean(axis=1)
    meanHA = HAspec.mean(axis=1)
    meanLA = LAspec.mean(axis=1)

    # define local function for calling snicar
    def simulate_albedo(rds, rho, dz, alg):
        params = c.namedtuple("params","rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, solzen")
        params.grain_rds = rds
        params.rho_layers = rho
        params.layer_type = [1,1]
        params.dz = dz
        params.mss_cnc_glacier_algae = alg
        params.solzen = 53
        albedo, BBA = call_snicar(params)
        return albedo, BBA


    # call snicar to generat esimulated spectrum for LA and HA using predetermined params

    SNICARalbedoLA, BBA = simulate_albedo([800,800],[800,800],[0.001,0.1],[0,0])
    SNICARalbedoHA, BBA = simulate_albedo([900,900],[800,800],[0.001,0.08],[0,0])
    SNICARalbedoLA = SNICARalbedoLA[15:230]
    SNICARalbedoHA = SNICARalbedoHA[15:230]


    # plot figure
    x = np.arange(350,2500,10)
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,5))
    ax1.plot(x,meanCI,linestyle='--',alpha = 0.4,label='Clean ice (mean)')
    ax1.plot(x,meanLA,linestyle='-.',alpha = 0.4,label='Algal ice (mean)')
    ax1.plot(x,SNICARalbedoLA,linestyle='dotted',alpha = 0.4,label='Clean ice (model)')
    ax1.fill_between(x,meanCI,SNICARalbedoLA,alpha=0.2)
    ax1.fill_between(x,SNICARalbedoLA,meanLA,color='k',alpha=0.2)
    ax1.set_xlim(350,1500), ax1.legend(loc='best')
    ax1.set_ylabel('Albedo'), ax1.set_xlabel('Wavelength (nm)')

    ax2.plot(x,meanCI,linestyle='--',alpha = 0.4,label='Clean ice (mean)')
    ax2.plot(x,meanHA,linestyle='-.',alpha = 0.4,label='Algal ice (mean)')
    ax2.plot(x,SNICARalbedoHA,linestyle='dotted',alpha = 0.4,label='Clean ice (model)')
    ax2.fill_between(x,meanCI,SNICARalbedoHA,alpha=0.2)
    ax2.fill_between(x,SNICARalbedoHA,meanHA,color='k',alpha=0.2)
    ax2.set_xlim(350,1500), ax2.legend(loc='best')
    ax2.set_ylabel('Albedo'), ax2.set_xlabel('Wavelength (nm)')
    fig.tight_layout()
    plt.savefig(str(savepath+'/BiovsPhysEffect.jpg'),dpi=300)
        
    # define incoming to calculate broadband albedo
    incoming = xr.open_dataset('/home/joe/Code/BioSNICAR_GO_PY/Data/Mie_files/480band/fsds/swnb_480bnd_sas_clr_SZA60.nc')
    incoming = incoming['flx_frc_sfc'].values
    incoming = incoming[15:230]

    # calculate broadband albedo of each case
    LA_BBA = np.sum(meanLA*incoming)/np.sum(incoming)
    HA_BBA = np.sum(meanHA*incoming)/np.sum(incoming)
    CI_BBA = np.sum(meanCI*incoming)/np.sum(incoming)
    CI2_BBA_LA = np.sum(SNICARalbedoLA*incoming)/np.sum(incoming)
    CI2_BBA_HA = np.sum(SNICARalbedoHA*incoming)/np.sum(incoming)

    # calculate change due to bio/phys as BBA difference
    delAbioLA = CI2_BBA_LA-LA_BBA
    delAphysLA = CI_BBA-CI2_BBA_LA
    delAbioHA = CI2_BBA_HA-HA_BBA
    delAphysHA = CI_BBA-CI2_BBA_HA


    return delAbioLA,delAphysLA,delAbioHA,delAphysHA



def build_LUT(solzen, dz, densities, radii, algae, wavelengths, save_LUT, apply_ARF, ARF_CI, ARF_HA, savepath):

    """
    generates LUTs used to invert BioSNICAR in RISA project

    params:
    ice_rds: fixed effective bubble radius for solid ice layers (default = 525)
    ice dens: fixed density for solid ice layers (default = 894)
    zeniths: range of solar zenith angles to loop over
    dz: thickness of each vertical layer
    densities: densities for top layer. Lower layers predicted by exponential model
    algae: mass mixing ratio of algae in top layer
    wavelengths: wavelength range, default is np.arange(0.2, 5, 0.01)
    save_LUT: Boolean to toggle saving to npy file
    savepath: directory to save LUT

    returns:
    WCthickLUT: for each index position in the spectraLUT, this holds the WC 
                thickness in the corresponding index position
    SpectraLUT: ND array containing 480element spectrum for each 
                dens/alg/zen combination

    return spectraLUT 


    """
    LUT = []

    @dask.delayed
    def run_sims(dens,rad,dz,alg,zen):

        params = c.namedtuple("params", "rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, solzen")
        params.rho_layers = [dens,dens]
        params.grain_rds = [rad,rad] # set equal to density
        params.layer_type = [1,1]
        params.dz = [0.001,dz]
        params.mss_cnc_glacier_algae = [alg,0]
        params.solzen = zen

        albedo, BBA = call_snicar(params)

        return albedo

    for z in np.arange(0,len(solzen),1):
        for i in np.arange(0,len(densities),1):
            for j in np.arange(0,len(radii),1):
                for p in np.arange(0,len(dz),1):
                    for q in np.arange(0,len(algae),1):
                        
                        albedo = run_sims(densities[i],radii[j],dz[p],algae[q],solzen[z])
                                
                        LUT.append(albedo)

    LUT = dask.compute(*LUT,scheduler='processes')
    LUT = np.array(LUT).reshape(len(solzen),len(densities),len(radii),len(dz),len(algae),len(wavelengths))

    # move the ARF application to new loop because dask compute objets are immutable
    # i.e. modifications to albedo must be done post-compute
    if apply_ARF:
        for z in np.arange(0,len(solzen),1):
            for i in np.arange(0,len(densities),1):
                for j in np.arange(0,len(radii),1):
                    for p in np.arange(0,len(dz),1):
                        for q in np.arange(0,len(algae),1):

                            if algae[q] > 5000:
                                
                                LUT[z,i,j,p,q,0:215] = LUT[z,i,j,p,q,0:215]*ARF_HA
                                
                            else:
                                LUT[z,i,j,p,q,0:215] = LUT[z,i,j,p,q,0:215]*ARF_CI


    if save_LUT:
        np.save(str(savepath+"LUT.npy"),LUT)

    return LUT



def inverse_model(field_data_fname,path_to_LUTs):

    spectra = pd.read_csv(field_data_fname)

    LUT_idx = [19, 26, 36, 40, 44, 48, 56, 131, 190]
    spectrum_idx = [140, 210, 315, 355, 390, 433, 515, 1260, 1840]

    algaeList = []
    grainList = []
    densityList =[]
    dzList = []
    filenames = []
    errorList = []
    zenithList = []


    for i in np.arange(0,len(spectra.columns),1):
        
        if i != 'wavelength':

            colname = spectra.columns[i]
            spectrum = np.array(spectra[colname])
            # calculate 2BDA index of field spectrum
            BDA2idx = np.array(spectrum)[360]/np.array(spectrum)[330]
            BDA2cells = abs(216000*BDA2idx-208600)

            print(BDA2cells)

            if (BDA2cells < 5000):
                LUT = np.load(str(path_to_LUTs+'LUT_1.npy'))
                dz = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2]
                densities = [600, 650, 700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [0, 2500, 5000]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT1")

            elif (BDA2cells >= 5000)&(BDA2cells < 10000):
                LUT = np.load(str(path_to_LUTs+'LUT_2.npy'))
                dz = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [7500, 10000, 12500]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT2")

            elif (BDA2cells >= 10000)&(BDA2cells < 15000):
                LUT = np.load(str(path_to_LUTs+'LUT_3.npy'))
                dz = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [10000, 12500, 15000, 17500]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT3")

            elif (BDA2cells >=15000)&(BDA2cells<20000):
                LUT = np.load(str(path_to_LUTs+'LUT_4.npy'))
                dz = [0.02, 0.03, 0.04, 0.05]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [15000, 20000, 25000]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT4")

            elif (BDA2cells >=20000)&(BDA2cells<25000):
                LUT = np.load(str(path_to_LUTs+'LUT_5.npy'))
                dz = [0.02, 0.03, 0.04, 0.05]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [20000, 22500, 25000, 27500, 30000]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT5")

            elif (BDA2cells >=25000)&(BDA2cells<30000):
                LUT = np.load(str(path_to_LUTs+'LUT_6.npy'))
                dz = [0.02, 0.03, 0.04, 0.05]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [25000, 27500, 30000, 32500, 35000]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT6")

            elif (BDA2cells >30000):
                LUT = np.load(str(path_to_LUTs+'LUT_7.npy'))
                dz = [0.02, 0.03, 0.04, 0.05]
                densities = [700, 750, 800, 850, 900]
                radii = [600, 700, 800, 900]
                algae = [35000, 37500, 40000, 45000]
                solzens = [40, 45, 50, 55]
                wavelengths = np.arange(0.2,5,0.01)
                print("LUT7")

            LUT = LUT[:,:,:,:,:,15:230]
            LUT = LUT.reshape(len(solzens)*len(densities)*len(radii)*len(dz)*len(algae),len(wavelengths[15:230]))
            
            # LUT = LUT[:,LUT_idx] # reduce wavelengths to only the 9 that match the S2 image
            # spectrum = spectrum[spectrum_idx]
            
            error_array = abs(LUT - spectrum[0:-1:10])

            mean_error = np.mean(error_array,axis=1)
            
            index = np.argmin(mean_error)

            param_idx = np.unravel_index(index,[len(solzens),len(densities),len(radii),len(dz),len(algae)])

            filenames.append(colname)
            zenithList.append(solzens[param_idx[0]])
            densityList.append(densities[param_idx[1]])
            grainList.append(radii[param_idx[2]])
            dzList.append(dz[param_idx[3]])
            algaeList.append(algae[param_idx[4]])
            errorList.append(np.min(mean_error))
            
    Out = pd.DataFrame(columns=['filename','density','grain','algae'])
    Out['filename'] = filenames
    Out['zenith'] = zenithList
    Out['density'] = densityList
    Out['grain'] = grainList
    Out['dz'] = dzList
    Out['algae'] = algaeList
    Out['spec_error'] = errorList

    return Out


def BDA2_of_field_samples():

    """
    2DBA index calculated from field samples after field spectra are averaged over S2 band 4 and 5 wavelengths
    weighted by the sensor spectral response function for each band. The index is then calculated as B5/B4
    and the cell concentration predicted using Wang et al's (2018) conversion equation.

    """

    spectra = pd.read_csv('/home/joe/Code/Remote_Ice_Surface_Analyser/Training_Data/HCRF_master_16171819.csv')
    
    # reformat LUT: flatten LUT from 3D to 2D array with one column per combination
    # of RT params, one row per wavelength
    
    responsefunc = pd.read_csv('/home/joe/Code/Remote_Ice_Surface_Analyser/S2SpectralResponse.csv')
    func04 = responsefunc['B4'].loc[(responsefunc['SR_WL']>650)&(responsefunc['SR_WL']<=680)]
    func05 = responsefunc['B5'].loc[(responsefunc['SR_WL']>698)&(responsefunc['SR_WL']<=713)]

    filenames = []
    Idx2DBAList =[]
    prd2DBAList = []
    Idx2DBA_S2List =[]
    prd2DBA_S2List = []
    Idx2DBA_Ideal_List = []
    prd2DBA_Ideal_List = []

    for i in np.arange(0,len(spectra.columns),1):
        
        if i != 'Wavelength':

            colname = spectra.columns[i]
            spectrum = np.array(spectra[colname])
            
            B04 = np.mean(spectrum[300:330] * func04)
            B05 = np.mean(spectrum[348:363] * func05)
            Idx2DBA = spectrum[355]/spectrum[315]
            prd2DBA = 10E-35 * Idx2DBA * np.exp(87.015*Idx2DBA)
            Idx2DBA_S2 = B05/B04
            prd2DBA_S2 = 10E-35 * Idx2DBA_S2 * np.exp(87.015*Idx2DBA_S2)
            Idx2DBA_Ideal = spectrum[360]/spectrum[330]
            prd2DBA_Ideal = 10E-35 * Idx2DBA_Ideal * np.exp(87.015*Idx2DBA_Ideal)
            
            filenames.append(colname)
            Idx2DBAList.append(Idx2DBA)
            prd2DBAList.append(prd2DBA)
            Idx2DBA_S2List.append(Idx2DBA_S2)
            prd2DBA_S2List.append(prd2DBA_S2)
            Idx2DBA_Ideal_List.append(Idx2DBA_Ideal)
            prd2DBA_Ideal_List.append(prd2DBA_Ideal)

    Out = pd.DataFrame()
    Out['filename'] = filenames
    Out['2DBAIdx'] = Idx2DBAList
    Out['2DBAPrediction'] = prd2DBAList
    Out['2DBA_S2Idx'] = Idx2DBA_S2List
    Out['2DBA_S2Prediction'] = prd2DBA_S2List
    Out['2DBAIdx_Ideal'] = Idx2DBA_Ideal_List
    Out['2DBAPrediction_Ideal'] = prd2DBA_Ideal_List

    return Out


def compare_predicted_and_measured(savepath, path_to_metadata):

    ## imports and data organisation
    import statsmodels.api as sm
    import numpy as np
    import pandas as pd


    DF = pd.read_csv(path_to_metadata)

    measured_cells = DF['measured_cells']
    modelled_cells = DF['algae_cells_inv_model_S2']
    BDA2_cells = DF['cells_BDA2_centre_wang']
    BDA2_idx = DF['BDA2_centre']

    
    ## regression models

    # Ordinary least squares regression
    model1 = sm.OLS(modelled_cells,measured_cells).fit()
    summary1 = model1.summary()
    test_x = [0,1000,5000,7500,10000,12500,15000,17500,20000,25000, 30000, 35000, 40000, 50000]
    ypred1 = model1.predict(test_x)

    # regress measured cells against band index 
    # use this to give predictive linear model 

    BDA2_PredModel = sm.OLS(measured_cells,sm.add_constant(BDA2_idx)).fit()
    BDA2_PredModel_r2 = np.round(BDA2_PredModel.rsquared,3)
    BDA2_PredModel_y = BDA2_PredModel.predict(sm.add_constant(BDA2_idx)) 
    
    # regress BDA2 predicted cells against measured cells
    model2 = sm.OLS(BDA2_PredModel_y,measured_cells).fit()
    summary2 = model2.summary()
    test_x = [0,1000,5000,7500,10000,12500,15000,17500,20000,25000, 30000, 35000, 40000, 50000]
    ypred2 = model2.predict(test_x)

    # multipanel figure
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(8,8))

    ax1.plot(measured_cells,color='k',marker = 'x',\
        linestyle='None',label='field-measured')
    ax1.plot(modelled_cells,color='b',marker = 'o', \
        markerfacecolor='None', alpha=0.6, linestyle = 'None',\
            label='RTM model prediction')
    ax1.plot(BDA2_PredModel_y,color='r', marker ='^', \
        markerfacecolor='r', alpha=0.3, linestyle = 'None',\
            label='new 2BDA model prediction')
    ax1.set_ylabel('Algal concentration (cells/mL)')    
    ax1.set_xticks(range(len(measured_cells)))
    ax1.set_xticklabels([])
    ax1.legend(loc='upper left')
    ax1.set_xlabel('Individual samples')
    ax1.set_ylim(0,65000)

    ax2.scatter(measured_cells, modelled_cells, marker='o',\
        facecolor ='None',color='b',alpha=0.6,\
        label='RTM\nr$^2$ = {}\np = {}'.format(np.round(model1.rsquared,3),\
            np.round(model1.pvalues[0],8)))
    ax2.plot(test_x, ypred1, linestyle='dotted',color='b',alpha=0.6)
    ax2.scatter(measured_cells, BDA2_PredModel_y, marker = '^',\
         facecolor='r',color='r', alpha=0.3, label='2BDA\nr$^2$ = {}\np = {}'.format(\
             np.round(model2.rsquared,3),\
        np.round(model2.pvalues[0],10)))
    ax2.plot(test_x, ypred2, linestyle = 'dashed', color='r',alpha=0.3)
    ax2.set_ylabel('Algal concentration,\n cells/mL (field)')
    ax2.set_xlabel('Algal concentration,\n clls/mL (predicted by model)')
    ax2.set_xlim(0,50000),ax2.set_ylim(0,60000)

    ax2.legend(loc='upper left')

    fig.tight_layout()

    savepath = '/home/joe/Code/Remote_Ice_Surface_Analyser/Manuscript/Figures'
    fig.savefig(str(savepath+'/measured_modelled_algae.png'),dpi=300)

    return
