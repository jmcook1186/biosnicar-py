"""
#####################################################################
################# BioSNICAR_GO DRIVER SCRIPT ########################

This script is used to configure the 2-stream radiative transfer
model BioSNICAR_GO. Here variable values are defined, the model called
and the results plotted.

NB. Setting Mie = 1, GO = 0 and algal impurities = 0 is equivalent to
running the original SNICAR model of Flanner et al. (2007, 2009)

NB: if using only granular layers, recommend using the faster Toon et al
tridiagonal matix solver (by setting TOON = True), however this will not
include any specular reflection components. If solid ice layers are
included in the ice/snow column, the ADDING-DOUBLING solver must be used
(i.e. ADD_DOUBLE = True).

Author: Joseph Cook, October 2020


TODO: add new laps according to 480-band optical property files (currently just one glacier alga
with dims 4 x 40)

TODO: Update README

TODO: create new algal optical properties that extend to 200nm and use full phenol spectrum
(requires generating new optical property files)

######################################################################
######################################################################

"""

from SNICAR_feeder import snicar_feeder
import matplotlib.pyplot as plt
import numpy as np

# set dir_base to the location of the BioSNICAR_GO_PY folder

dir_base = '/home/joe/Code/BioSNICAR_GO_PY/'

##############################
# 1) Choose plot/print options
###############################


show_figs = True # toggle to display spectral albedo figure
save_figs = True # toggle to save spectral albedo figure to file
savepath = dir_base # base path for saving figures
print_BBA = True # toggle to print broadband albedo to terminal
print_band_ratios = False # toggle to print various band ratios to terminal
smooth = True # apply optional smoothing function (Savitzky-Golay filter)
window_size = 11 # if applying smoothing filter, define window size
poly_order = 3 # if applying smoothing filter, define order of polynomial

##################################################################
# 2) CHOOSE METHOD FOR DETERMINING OPTICAL PROPERTIES OF ICE GRAINS
# for small spheres choose Mie, for hexagonal plates or columns of any size
# choose GeometricOptics
#####################################################################

MIE = True
GO = False

# If Mie = True, select solver (only Toon available for GO mode).
TOON = False # toggle Toon et al tridiagonal matrix solver
ADD_DOUBLE = True # toggle adding-doubling solver


######################################
## 3. RADIATIVE TRANSFER CONFIGURATION
#######################################

DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
solzen   = 20      # if DIRECT give solar zenith angle (degrees from 0 = nadir, 90 = horizon)
rf_ice = 2        # define source of ice refractive index data. 0 = Warren 1984, 1 = Warren 2008, 2 = Picard 2016

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
incoming_i = 0


#############################################
## 4. SET PHYSICAL PROPERTIES OF THE ICE/SNOW
#############################################

# grain shapes are only active in Mie scattering mode - GO assumes hexagonal columns.
# snw_shp can be 0 = sphere, 1 = spheroid, 2 = hexagonal plate, 3= koch snowflake
# shp_fctr = ratio of nonspherical grain effective radii to that of equal-volume sphere
    
    # 0=use recommended default value (He et al. 2017);
    # use user-specified value (between 0 and 1)
    # only activated when sno_shp > 1 (i.e. nonspherical)

dz = [0.001, 0.01, 1, 1, 10] # thickness of each vertical layer (unit = m)
nbr_lyr = len(dz)  # number of snow layers
layer_type = [0,0,1,1,1]
R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
rho_snw = [600, 600, 894, 894, 894] # density of each layer (unit = kg m-3)
rds_snw = [1500,1500,550,550,550] # effective grain radius of snow/bubbly ice
rwater = [0, 0, 0, 0, 0] # if  using Mie calculations, add radius of optional liquid water coating
snw_shp =[0,0,0,0,0] # grain shape(He et al. 2016, 2017)
shp_fctr = [0,0,0,0,0] # shape factor (ratio of aspherical grain radii to that of equal-volume sphere)
snw_ar = [0,0,0,0,0] # aspect ratio (ratio of width to length)

# if using GeometricOptics, set side_length and depth
side_length = [10000,10000,10000,10000,10000] 
depth = [20000,20000,20000,20000,20000]

#######################################
## 5) SET LAP CHARACTERISTICS
#######################################

nbr_aer = 25 # Define total number of different LAPs/aerosols in model

# set filename stubs
stb1 = 'RealPhenol_algae_geom_' # %name stub 1
stb2 = '.nc'  # file extension
wrkdir2 = str(dir_base + '/Data/Algal_Optical_Props/') # working directory
snw_stb1 = 'snw_alg_' # name stub for snow algae

# CHOOSE DIMENSIONS OF GLACIER ALGAE 1
algae_r = 6 # algae radius
algae_l = 60 # algae length
#glacier_algae1 = str(wrkdir2+stb1+str(algae_r)+'_'+str(algae_l)+stb2) # create filename string
glacier_algae1 = str(wrkdir2+'RealPhenol_algae_geom_{}_{}.nc'.format(algae_r,algae_l))

# CHOOSE DIMENSIONS OF GLACIER ALGAE 2
algae2_r = 2 # algae radius
algae2_l = 10 # algae length
glacier_algae2 = str(wrkdir2+stb1+str(algae2_r)+'_'+str(algae2_l)+stb2) # create filename string

# SET UP IMPURITY MIXING RATIOS
# PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)

for x in [100000]:
    
    mss_cnc_soot1 = [0,0,0,0,0]    # uncoated black carbon (Bohren and Huffman, 1983)
    mss_cnc_soot2 = [0,0,0,0,0]    # coated black carbon (Bohren and Huffman, 1983)
    mss_cnc_brwnC1 = [0,0,0,0,0]   # uncoated brown carbon (Kirchstetter et al. (2004).)
    mss_cnc_brwnC2 = [0,0,0,0,0]   # sulfate-coated brown carbon (Kirchstetter et al. (2004).)
    mss_cnc_dust1 = [0,0,0,0,0]    # dust size 1 (r=0.05-0.5um) (Balkanski et al 2007)
    mss_cnc_dust2 = [0,0,0,0,0]    # dust size 2 (r=0.5-1.25um) (Balkanski et al 2007)
    mss_cnc_dust3 = [0,0,0,0,0]    # dust size 3 (r=1.25-2.5um) (Balkanski et al 2007)
    mss_cnc_dust4 = [0,0,0,0,0]    # dust size 4 (r=2.5-5.0um)  (Balkanski et al 2007)
    mss_cnc_dust5 = [0,0,0,0,0]    # dust size 5 (r=5.0-50um)  (Balkanski et al 2007)
    mss_cnc_ash1 = [0,0,0,0,0]    # volcanic ash size 1 (r=0.05-0.5um) (Flanner et al 2014)
    mss_cnc_ash2 = [0,0,0,0,0]    # volcanic ash size 2 (r=0.5-1.25um) (Flanner et al 2014)
    mss_cnc_ash3 = [0,0,0,0,0]    # volcanic ash size 3 (r=1.25-2.5um) (Flanner et al 2014)
    mss_cnc_ash4 = [0,0,0,0,0]    # volcanic ash size 4 (r=2.5-5.0um) (Flanner et al 2014)
    mss_cnc_ash5 = [0,0,0,0,0]    # volcanic ash size 5 (r=5.0-50um) (Flanner et al 2014)
    mss_cnc_Skiles_dust1 = [0,0,0,0,0]    # Colorado dust size 1 (Skiles et al 2017)
    mss_cnc_Skiles_dust2 = [0,0,0,0,0]    # Colorado dust size 2 (Skiles et al 2017)
    mss_cnc_Skiles_dust3 = [0,0,0,0,0]    # Colorado dust size 3 (Skiles et al 2017)
    mss_cnc_Skiles_dust4 = [0,0,0,0,0]  # Colorado dust size 4 (Skiles et al 2017)
    mss_cnc_Skiles_dust5 = [0,0,0,0,0]  # Colorado dust size 5 (Skiles et al 2017)
    mss_cnc_GreenlandCentral1 = [0,0,0,0,0] # Greenland Central dust size 1 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral2 = [0,0,0,0,0] # Greenland Central dust size 2 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral3 = [0,0,0,0,0] # Greenland Central dust size 3 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral4 = [0,0,0,0,0] # Greenland Central dust size 4 (Polashenski et al 2015)
    mss_cnc_GreenlandCentral5 = [0,0,0,0,0] # Greenland Central dust size 5 (Polashenski et al 2015)
    mss_cnc_snw_alg = [0,0,0,0,0]    # Snow Algae (spherical, C nivalis) (Cook et al. 2017)
    mss_cnc_glacier_algae = [0,0,0,0,0]    # glacier algae type1 (Cook et al. 2020)


    ##########################################################################
    ################## CALL FUNCTIONS AND PLOT OUTPUTS #######################
    ##########################################################################

    # SET FILE NAMES CONTAINING OPTICAL PARAMETERS FOR ALL IMPURITIES:

    FILE_soot1  = '/480band/lap/mie_sot_ChC90_dns_1317.nc'
    FILE_soot2  = '/480band/lap/miecot_slfsot_ChC90_dns_1317.nc'
    FILE_brwnC1 = '/480band/lap/brC_Kirch_BCsd.nc'
    FILE_brwnC2 = '/480band/lap/brC_Kirch_BCsd_slfcot.nc'
    FILE_dust1  = '/480band/lap/dust_balkanski_central_size1.nc'
    FILE_dust2  = '/480band/lap/dust_balkanski_central_size2.nc'
    FILE_dust3  = '/480band/lap/dust_balkanski_central_size3.nc'
    FILE_dust4  = '/480band/lap/dust_balkanski_central_size4.nc'
    FILE_ash1  = '/480band/lap/volc_ash_eyja_central_size1.nc'
    FILE_ash2 = '/480band/lap/volc_ash_eyja_central_size2.nc'
    FILE_ash3 = '/480band/lap/volc_ash_eyja_central_size3.nc'
    FILE_ash4 = '/480band/lap/volc_ash_eyja_central_size4.nc'
    FILE_ash5 = '/480band/lap/volc_ash_eyja_central_size5.nc'
    FILE_Skiles_dust1 = '/480band/lap/dust_skiles_size1.nc'
    FILE_Skiles_dust2 = '/480band/lap/dust_skiles_size2.nc'
    FILE_Skiles_dust3 = '/480band/lap/dust_skiles_size3.nc'
    FILE_Skiles_dust4 = '/480band/lap/dust_skiles_size4.nc'
    FILE_Skiles_dust5 = '/480band/lap/dust_skiles_size5.nc'
    FILE_GreenlandCentral1 = '/480band/lap/dust_greenland_central_size1.nc'
    FILE_GreenlandCentral2 = '/480band/lap/dust_greenland_central_size2.nc'
    FILE_GreenlandCentral3 = '/480band/lap/dust_greenland_central_size3.nc'
    FILE_GreenlandCentral4 = '/480band/lap/dust_greenland_central_size4.nc'
    FILE_GreenlandCentral5  = '/480band/lap/dust_greenland_central_size5.nc'
    FILE_snw_alg  = '/480band/lap/snw_alg_r025um_chla020_chlb025_cara150_carb140.nc'
    FILE_glacier_algae = '/480band/lap/Glacier_Algae_480.nc'


    #########################################################
    # Error catching: invalid combinations of input variables
    #########################################################

    if MIE == True and GO == True:

        raise ValueError("ERROR: BOTH MIE AND GO MODES SELECTED: PLEASE CHOOSE ONE")

    elif TOON == True and ADD_DOUBLE == True:

        raise ValueError("ERROR: BOTH SOLVERS SELECTED: PLEASE CHOOSE EITHER TOON OR ADD_DOUBLE")

    elif TOON == True and solzen < 40:
        
        raise ValueError("INVALID SOLAR ANGLE: solzen outside valid range for Toon solver - ue AD mode or solzen >= 40.")
    
    elif np.sum(layer_type) < 1 and ADD_DOUBLE==True:
        # just warn user but let program continue - in some cases 
        # AD method preferable (stable over complete range of SZA)
        print("\nWARNING:")
        print("There are no solid ice layers included - you might prefer the faster matrix inversion method.")
        print("Toggle TOON=True and ADD_DOUBLE=False to use it.\n")

    elif np.sum(layer_type) > 0 and TOON == True:

        raise ValueError("There are ice layers in the model - please use the adding-doubling solver")

    if np.sum(mss_cnc_snw_alg) != 0:
        # remind user that snow algae optical properties have not yet been empirically validated
        print("WARNING: you are using snow algae as an impurity in the model.")
        print("the optical properties for these algae are theoretical.") 
        print("They were constructed from literature values for pigmentation, refractive indices and cell size")
        print("They have not yet been validated empirically.")

    #######################################
    # IF NO INPUT ERRORS --> FUNCTION CALLS
    #######################################

        
    [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt] =\
    snicar_feeder(MIE, GO, dir_base, rf_ice, incoming_i, DIRECT, layer_type,\
    APRX_TYP, DELTA, solzen, TOON, ADD_DOUBLE, R_sfc, dz, rho_snw, rds_snw,\
    side_length, depth, rwater, nbr_lyr, nbr_aer, snw_shp, shp_fctr, snw_ar,\
    mss_cnc_soot1, mss_cnc_soot2, mss_cnc_brwnC1, mss_cnc_brwnC2, mss_cnc_dust1,\
    mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_ash2,\
    mss_cnc_ash3, mss_cnc_ash4, mss_cnc_ash5, mss_cnc_Skiles_dust1, mss_cnc_Skiles_dust2,\
    mss_cnc_Skiles_dust3, mss_cnc_Skiles_dust4, mss_cnc_Skiles_dust5, mss_cnc_GreenlandCentral1,\
    mss_cnc_GreenlandCentral2, mss_cnc_GreenlandCentral3, mss_cnc_GreenlandCentral4,\
    mss_cnc_GreenlandCentral5, mss_cnc_snw_alg, mss_cnc_glacier_algae, FILE_soot1,\
    FILE_soot2, FILE_brwnC1, FILE_brwnC2, FILE_dust1, FILE_dust2, FILE_dust3, FILE_dust4,\
    FILE_ash1, FILE_ash2, FILE_ash3, FILE_ash4, FILE_ash5, FILE_Skiles_dust1, FILE_Skiles_dust2,\
    FILE_Skiles_dust3, FILE_Skiles_dust4, FILE_Skiles_dust5, FILE_GreenlandCentral1,\
    FILE_GreenlandCentral2, FILE_GreenlandCentral3, FILE_GreenlandCentral4, FILE_GreenlandCentral5,\
    FILE_snw_alg, FILE_glacier_algae)

   
    ########################
    # PLOTTING AND PRINTING
    ########################

    if smooth:
        from scipy.signal import savgol_filter
        yhat = savgol_filter(albedo, window_size, poly_order)
        albedo = yhat

    if print_band_ratios:

        I2DBA = albedo[40]/albedo[36]
        I3DBA = (albedo[36] - albedo[40]) / albedo[45]
        NDCI = ((albedo[40]-albedo[38])-(albedo[45]-albedo[38]))*((albedo[40]-albedo[38])/(albedo[45]-albedo[38]))
        MCI = (albedo[40]-albedo[36])/(albedo[40]+albedo[36])
        II = np.log(albedo[26])/np.log(albedo[56])

        print("\nINDEX VALUES")
        print("2DBA Index: ",I2DBA)
        print("3DBA index: ", I3DBA)
        print("NDCI index: ", NDCI)
        print("MCI index: ", MCI)
        print("Impurity Index: ", II)

    if print_BBA:

        print('\nBROADBAND ALBEDO = ', BBA)

    #### PLOT ALBEDO ######
    plt.plot(wvl, albedo)

plt.ylabel('ALBEDO'), plt.xlabel('WAVELENGTH (microns)'), plt.xlim(0.2,1.8),
plt.ylim(0,1), plt.axvline(x = 0.68,color='g',linestyle='dashed')

print(len(albedo))

if save_figs:
    plt.savefig(str(savepath+"spectral_albedo.png"))
    plt.show()

if show_figs:
    plt.show()