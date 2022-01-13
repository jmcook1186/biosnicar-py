"""
Joseph Cook, Jan 2021:

Functions for building lookup table and extracting band ratio values for 
snicar model of the weathering crust. Separate from other functions in RISA
projct because these functions call to SNICAR and to SNICAR parameterisation
scripts and therefore need to be run in the BioSNICAR_GO_PY folder and
BioSNICAR_py conda environment.

"""
from snicar_feeder import snicar_feeder
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import collections


def build_LUT(ice_rds,\
    ice_dens,\
    zeniths,\
    dz,\
    densities,algae,wavelengths, save_LUT, savepath):

    """
    generates LUT used to invert BioSNICAR in RISA project

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


    NOTE: r_eff is not included in the LUT. The forward modelling showed that
    the r_eff was usually equal to the density value and when it was not equal 
    it was +/- 50. Setting the r_eff == density increased the mean absolute error
    in our field comparisons by ~0.001. This heuristic reduces the d.o.f of the LUT
    and allows more options for density, algae and thickness to be included.

    """

    # LUT dims = [density, r_eff, thickness, algae, wavelengths]
    spectraLUT = np.zeros(shape = (len(densities), 3, len(dz),len(algae),len(wavelengths)))

    for i in np.arange(0,len(zeniths),1):
        for j in np.arange(0,len(densities),1):
            for p in np.arange(0,len(dz),1):
                for q in np.arange(0,len(algae),1):
                    
                   
                    params = collections.namedtuple("params","rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, solzen")
                    params.grain_rds = [densities[j],densities[j]] # set equal to density
                    params.rho_layers = [densities[j],densities[j]]
                    params.layer_type = [1,1]
                    params.dz = [0.001,dz[p]]
                    params.mss_cnc_glacier_algae = [algae[q],0]
                    params.solzen = 60

                    albedo, BBA = call_snicar(params)

                    spectraLUT[j,p,q,:] = albedo

        if save_LUT:
            zen_name = int(np.round((np.cos(zeniths[i] * (np.pi / 180))),2)*100)
            np.save(str(savepath+"Spec_LUT_{}.npy".format(zen_name)),spectraLUT)
            np.save(str(savepath+"WC_LUT_{}.npy".format(zen_name)),WCthickLUT)

    return spectraLUT 



def call_snicar(params):
    
# set dir_base to the location of the BioSNICAR_GO_PY folder

    dir_base = '/home/joe/Code/BioSNICAR_GO_PY/'
    savepath = dir_base # base path for saving figures
       
    TOON = False # toggle Toon et al tridiagonal matrix solver
    ADD_DOUBLE = True # toggle adding-doubling solver
    layer_type = params.layer_type
    DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
    solzen   = params.solzen      # if DIRECT give solar zenith angle (degrees from 0 = nadir, 90 = horizon)
    rf_ice = 2        # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016
    incoming_i = 4
    nbr_lyr = len(params.dz)  # number of snow layers
    R_sfc = 0.1 # reflectance of underlying surface - set across all wavelengths
    rwater = [0]*len(params.dz) # if  using Mie calculations, add radius of optional liquid water coating
    grain_shp =[0]*len(params.dz) # grain shape(He et al. 2016, 2017)
    shp_fctr = [0]*len(params.dz) # shape factor (ratio of aspherical grain radii to that of equal-volume sphere)
    grain_ar = [0]*len(params.dz) # aspect ratio (ratio of width to length)
    side_length = 0
    depth=0
    grain_rds = params.grain_rds
    rho_layers = params.rho_layers
    dz = params.dz
    
    mss_cnc_soot1 = [0]*len(params.dz)    # uncoated black carbon (Bohren and Huffman, 1983)
    mss_cnc_soot2 = [0]*len(params.dz)    # coated black carbon (Bohren and Huffman, 1983)
    mss_cnc_brwnC1 = [0]*len(params.dz)   # uncoated brown carbon (Kirchstetter et al. (2004).)
    mss_cnc_brwnC2 = [0]*len(params.dz)   # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    mss_cnc_dust1 = [0]*len(params.dz)    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    mss_cnc_dust2 = [0]*len(params.dz)    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    mss_cnc_dust3 = [0]*len(params.dz)    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    mss_cnc_dust4 = [0]*len(params.dz)    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    mss_cnc_dust5 = [0]*len(params.dz)    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    mss_cnc_ash1 = [0]*len(params.dz)    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    mss_cnc_ash2 = [0]*len(params.dz)    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    mss_cnc_ash3 = [0]*len(params.dz)    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    mss_cnc_ash4 = [0]*len(params.dz)    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    mss_cnc_ash5 = [0]*len(params.dz)    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    mss_cnc_ash_st_helens = [0]*len(params.dz)   # ash from Mount Saint Helen's
    mss_cnc_Skiles_dust1 = [0]*len(params.dz)    # Colorado dust size 1 (Skiles et al 2017)
    mss_cnc_Skiles_dust2 = [0]*len(params.dz)    # Colorado dust size 2 (Skiles et al 2017)
    mss_cnc_Skiles_dust3 = [0]*len(params.dz)    # Colorado dust size 3 (Skiles et al 2017)
    mss_cnc_Skiles_dust4 = [0]*len(params.dz)  # Colorado dust size 4 (Skiles et al 2017)
    mss_cnc_Skiles_dust5 = [0]*len(params.dz)  # Colorado dust size 5 (Skiles et al 2017)
    mss_cnc_GreenlandCentral1 = [0]*len(params.dz) # Greenland Central dust size 1 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral2 = [0]*len(params.dz) # Greenland Central dust size 2 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral3 = [0]*len(params.dz) # Greenland Central dust size 3 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral4 = [0]*len(params.dz) # Greenland Central dust size 4 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral5 = [0]*len(params.dz) # Greenland Central dust size 5 (Polashenski et al 2015)
    mss_cnc_Cook_Greenland_dust_L = [0]*len(params.dz)
    mss_cnc_Cook_Greenland_dust_C = [0]*len(params.dz)
    mss_cnc_Cook_Greenland_dust_H = [0]*len(params.dz)
    mss_cnc_snw_alg = [0]*len(params.dz)    # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
    mss_cnc_glacier_algae = params.mss_cnc_glacier_algae   # glacier algae type1 (Cook et al. 2020)

    nbr_aer = 30

    # Set names of files containing the optical properties of these LAPs:
    FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
    FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
    FILE_brwnC1 = 'brC_Kirch_BCsd.nc'
    FILE_brwnC2 = 'brC_Kirch_BCsd_slfcot.nc'
    FILE_dust1  = 'dust_balkanski_central_size1.nc'
    FILE_dust2  = 'dust_balkanski_central_size2.nc'
    FILE_dust3  = 'dust_balkanski_central_size3.nc'
    FILE_dust4  = 'dust_balkanski_central_size4.nc'
    FILE_dust5 = 'dust_balkanski_central_size5.nc'
    FILE_ash1  = 'volc_ash_eyja_central_size1.nc'
    FILE_ash2 = 'volc_ash_eyja_central_size2.nc'
    FILE_ash3 = 'volc_ash_eyja_central_size3.nc'
    FILE_ash4 = 'volc_ash_eyja_central_size4.nc'
    FILE_ash5 = 'volc_ash_eyja_central_size5.nc'
    FILE_ash_st_helens = 'volc_ash_mtsthelens_20081011.nc'
    FILE_Skiles_dust1 = 'dust_skiles_size1.nc'
    FILE_Skiles_dust2 = 'dust_skiles_size2.nc'
    FILE_Skiles_dust3 = 'dust_skiles_size3.nc'
    FILE_Skiles_dust4 = 'dust_skiles_size4.nc'
    FILE_Skiles_dust5 = 'dust_skiles_size5.nc'
    FILE_GreenlandCentral1 = 'dust_greenland_central_size1.nc'
    FILE_GreenlandCentral2 = 'dust_greenland_central_size2.nc'
    FILE_GreenlandCentral3 = 'dust_greenland_central_size3.nc'
    FILE_GreenlandCentral4 = 'dust_greenland_central_size4.nc'
    FILE_GreenlandCentral5  = 'dust_greenland_central_size5.nc'
    FILE_Cook_Greenland_dust_L = 'dust_greenland_Cook_LOW_20190911.nc'
    FILE_Cook_Greenland_dust_C = 'dust_greenland_Cook_CENTRAL_20190911.nc'
    FILE_Cook_Greenland_dust_H = 'dust_greenland_Cook_HIGH_20190911.nc'
    FILE_snw_alg  = 'snw_alg_r025um_chla020_chlb025_cara150_carb140.nc'
    FILE_glacier_algae = 'Glacier_Algae_480.nc'


    #######################################
    # IF NO INPUT ERRORS --> FUNCTION CALLS
    #######################################

        
    [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt] =\
    snicar_feeder(dir_base,\
    rf_ice, incoming_i, DIRECT, layer_type,\
    APRX_TYP, DELTA, solzen, TOON, ADD_DOUBLE, R_sfc, dz, rho_layers, grain_rds,\
    side_length, depth, rwater, nbr_lyr, nbr_aer, grain_shp, shp_fctr, grain_ar,\
    mss_cnc_soot1, mss_cnc_soot2, mss_cnc_brwnC1, mss_cnc_brwnC2, mss_cnc_dust1,\
    mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_dust5, mss_cnc_ash1, mss_cnc_ash2,\
    mss_cnc_ash3, mss_cnc_ash4, mss_cnc_ash5, mss_cnc_ash_st_helens, mss_cnc_Skiles_dust1, mss_cnc_Skiles_dust2,\
    mss_cnc_Skiles_dust3, mss_cnc_Skiles_dust4, mss_cnc_Skiles_dust5, mss_cnc_GreenlandCentral1,\
    mss_cnc_GreenlandCentral2, mss_cnc_GreenlandCentral3, mss_cnc_GreenlandCentral4,\
    mss_cnc_GreenlandCentral5, mss_cnc_Cook_Greenland_dust_L, mss_cnc_Cook_Greenland_dust_C,\
    mss_cnc_Cook_Greenland_dust_H, mss_cnc_snw_alg, mss_cnc_glacier_algae, FILE_soot1,\
    FILE_soot2, FILE_brwnC1, FILE_brwnC2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4, FILE_dust5,\
    FILE_ash1, FILE_ash2, FILE_ash3, FILE_ash4, FILE_ash5, FILE_ash_st_helens, FILE_Skiles_dust1, FILE_Skiles_dust2,\
    FILE_Skiles_dust3, FILE_Skiles_dust4, FILE_Skiles_dust5, FILE_GreenlandCentral1,\
    FILE_GreenlandCentral2, FILE_GreenlandCentral3, FILE_GreenlandCentral4, FILE_GreenlandCentral5,\
    FILE_Cook_Greenland_dust_L, FILE_Cook_Greenland_dust_C, FILE_Cook_Greenland_dust_H, FILE_snw_alg, FILE_glacier_algae)

    return albedo, BBA


##################################
# SET VALUES
##################################

savepath = "/home/joe/Code/BioSNICAR_GO_PY/"

ice_rds = 850
ice_dens = 900
dz = [0.02, 0.04]#[0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3]
densities = [500, 800]#[500, 550, 600, 650, 700, 750, 800, 850, 900, 910]
algae = [0,10000]#[0, 2500, 5000, 7500, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000]:
zeniths = [37]#[37, 45, 53, 60] # = coszen 80, 70, 60, 50

wavelengths = np.arange(0.2,5,0.01)
save_LUT= True


####################################
# CALL FUNCS
####################################

spectraLUT = build_LUT(ice_rds,ice_dens,zeniths,dz,densities,algae,wavelengths,save_LUT,savepath)

