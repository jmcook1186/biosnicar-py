"""
#####################################################################
################# BioSNICAR_GO DRIVER SCRIPT ########################

This script is used to configure the 2-stream radiative transfer
model BioSNICAR_GO. Here variable values are defined, the model called
and the results plotted. 

NB. Setting Mie = 1, GO = 0 and algal impurities = 0 is equivalent to
running the original SNICAR model of Flanner et al. (2007, 2009)

Author: Joseph Cook, October 2019

######################################################################
######################################################################

"""
from snicar8d_mie import snicar8d_mie
from snicar8d_GO import snicar8d_GO
import matplotlib.pyplot as plt
import numpy as np

##############################
# 1) Choose plot/print options

show_figs = True # toggle to display spectral albedo figure
save_figs = True # toggle to save spectral albedo figure to file
savepath = "/home/joe/Code/BioSNICAR_GO_PY/" # base path for saving figures
print_BBA = True # toggle to print broadband albedo to terminal
print_band_ratios = False # toggle to print various band ratios to terminal
smooth = True # apply optional smoothing function (Savitzky-Golay filter)
window_size = 11 # if applying smoothing filter, define window size
poly_order = 3 # if applying smoothing filter, define order of polynomial

##################################################################
# 2) CHOOSE METHOD FOR DETERMINING OPTICAL PROPERTIES OF ICE GRAINS
# for small spheres choose Mie, for hexagonal plates or columns of any size
# choose GeometricOptics

Mie = False
GeometricOptics = True


######################################
## 3. RADIATIVE TRANSFER CONFIGURATION

DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 3        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
coszen   = 0.3      # if DIRECT give cosine of solar zenith angle 

##################################################################
# TOGGLE SPECULAR REFLECTION ON/OFF and DEFINE ROUGHNESS

specular_reflection = True # apply specular reflection at upper surface
smoothness = 1 # degree of smoothness between 0-1 - scales fresnel coefficient

#############################################
## 4. SET PHYSICAL PROPERTIES OF THE ICE/SNOW

dz = [0.02, 0.02, 0.02, 0.02, 0.2] # thickness of each vertical layer (unit = m)
nbr_lyr = len(dz)  # number of snow layers
R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
rho_snw = [700, 700, 700, 700, 700] # density of each layer (unit = kg m-3)

#SET ICE GRAIN DIMENSIONS
# if using Mie optical properties, set rds_snw and 
# an optional coated liquid water sphere
rds_snw = [2000,2000,2000,2000,2000]
rwater = [0, 0, 0, 0, 0]

# if using GeometricOptics, set side_length and depth
side_length = [15000,15000,15000,15000,15000] 
depth = [15000,15000,15000,15000,15000]

##############################################
## 5) SET LAP CHARACTERISTICS

nbr_aer = 16 # Define total number of different LAPs/aerosols in model

# set filename stubs
stb1 = 'algae_geom_' # %name stub 1
stb2 = '.nc'  # file extension
wrkdir2 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/' # working directory
snw_stb1 = 'snw_alg_' # name stub for snow algae

# CHOOSE DIMENSIONS OF GLACIER ALGAE 1
algae_r = 6 # algae radius
algae_l = 60 # algae length
#glacier_algae1 = str(wrkdir2+stb1+str(algae_r)+'_'+str(algae_l)+stb2) # create filename string
glacier_algae1 = str(wrkdir2+'RealPhenol_algae_geom_6_60.nc')

# CHOOSE DIMENSIONS OF GLACIER ALGAE 2
algae2_r = 6 # algae radius
algae2_l = 20 # algae length
glacier_algae2 = str(wrkdir2+stb1+str(algae2_r)+'_'+str(algae2_l)+stb2) # create filename string

# CHOOSE SNOW ALGAE DIAMETER
snw_algae_r = 1 # snow algae diameter
snw_alg = str(wrkdir2+snw_stb1+str(snw_algae_r)+stb2) # create filename string

# SET UP IMPURITY MIXING RATIOS
# PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)

for x in [2000]:
    
    mss_cnc_soot1 = [0,0,0,0,0]    # uncoated black carbon
    mss_cnc_soot2 = [0,0,0,0,0]    # coated black carbon
    mss_cnc_dust1 = [0,0,0,0,0]    # global average dust 1
    mss_cnc_dust2 = [0,0,0,0,0]    # global average dust 2
    mss_cnc_dust3 = [0,0,0,0,0]    # global average dust 3
    mss_cnc_dust4 = [0,0,0,0,0]    # global average dust 4
    mss_cnc_ash1 = [0,0,0,0,0]    # volcanic ash species 1
    mss_cnc_GRISdust1 = [0,0,0,0,0]    # GRIS dust 1 (Cook et al. 2019 "mean")
    mss_cnc_GRISdust2 = [0,0,0,0,0]    # GRIS dust 2 (Cook et al. 2019 HIGH)
    mss_cnc_GRISdust3 = [0,0,0,0,0]    # GRIS dust 3 (Cook et al. 2019 LOW)
    mss_cnc_GRISdustP1 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: low hematite)
    mss_cnc_GRISdustP2 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
    mss_cnc_GRISdustP3 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
    mss_cnc_snw_alg = [0,0,0,0,0]    # Snow Algae (spherical, C nivalis)
    mss_cnc_glacier_algae1 = [0,0,0,0,0]    # glacier algae type1
    mss_cnc_glacier_algae2 = [0,0,0,0,0]    # glacier algae type2

    ##########################################################################
    ################## CALL FUNCTIONS AND PLOT OUTPUTS #######################

    # SET FILE NAMES CONTAINING OPTICAL PARAMETERS FOR ALL IMPURITIES:

    FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
    FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
    FILE_dust1  = 'aer_dst_bln_20060904_01.nc'
    FILE_dust2  = 'aer_dst_bln_20060904_02.nc'
    FILE_dust3  = 'aer_dst_bln_20060904_03.nc'
    FILE_dust4  = 'aer_dst_bln_20060904_04.nc'
    FILE_ash1  = 'volc_ash_mtsthelens_20081011.nc'
    FILE_GRISdust1 = 'dust_greenland_Cook_CENTRAL_20190911.nc'
    FILE_GRISdust2 = 'dust_greenland_Cook_HIGH_20190911.nc'
    FILE_GRISdust3 = 'dust_greenland_Cook_LOW_20190911.nc'
    FILE_GRISdustP1 = 'dust_greenland_L_20150308.nc'
    FILE_GRISdustP2 = 'dust_greenland_C_20150308.nc'
    FILE_GRISdustP3 = 'dust_greenland_H_20150308.nc'
    FILE_snw_alg  = snw_alg # snow algae (c nivalis)
    FILE_glacier_algae1 = glacier_algae1 # Glacier algae
    FILE_glacier_algae2 = glacier_algae2 # Glacier algae

    # Error catching: invalid combinations of inut variables

    if (DIRECT != 1) & (specular_reflection == True):
        raise("ERROR: Specular reflection applied to diffuse incident radiation \n Please select\
        direct beam and set beam angle or toggle specular reflection OFF.")
    
    if (specular_reflection == True) & ((coszen <0.3)or(coszen>0.9)):
        raise("ERROR: outside valid range of solar zenith angle")

    if Mie == True and GeometricOptics == True:

        print("ERROR: BOTH MIE AND GO MODES SELECTED: PLEASE CHOOSE ONE")

    elif Mie == True:

        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_mie(DIRECT, specular_reflection, smoothness, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, rds_snw, rwater, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

    elif GeometricOptics == True:

        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(DIRECT, specular_reflection, smoothness, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

    else:
        
        print("NEITHER MIE NOR GO MODE SELECTED: PLEASE CHOOSE ONE")
    
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
        print('BROADBAND ALBEDO = ', BBA)

    #### PLOT ALBEDO ######
    plt.plot(wvl, albedo)

plt.ylabel('ALBEDO'), plt.xlabel('WAVELENGTH (microns)'), plt.xlim(0.3,1.5),
plt.ylim(0,1), plt.axvline(x = 0.68,color='g',linestyle='dashed')

if save_figs:
    plt.savefig(str(savepath+"spectral_albedo.png"))
    plt.show()

if show_figs:
    plt.show()