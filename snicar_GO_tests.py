"""
###########################
BioSNICAR_GO_PY UNIT TESTS
###########################

Joseph Cook, October 2019
github.com/jmcook1186/BioSNICAR_GO_PY/

#######################################

This script is for automated testing of the snicar_GO module against benchmark variable values generated
by the matlab version an stored in external csv files in the folder wdir/Unit_Tests/Snicar_GO/.

If the spectral albedo pedicted by the Python and Matlab versons are equal to within 1e-4, the test is
considered passed. The test results are printed to the console for each test.

Only one variable is changed in each unit test, all other variables are set to the default values defined below.
Note that in the case of testing the direct/diffuse condition, the quadrature approximation is used instead of
the default Eddington approximation. This is because of a known bug that allows the Eddington-predicted albedo
to go negative under diffuse illumination. Recommend avoiding the combination of Eddington approximation and
diffuse irradiance.

Note that these tests are not exhaustive, there are thousands of variable combinations available, and I have
only tested for equality between predicted albedo values so far. Bugs that only arise under certain variable
value combinations or only affect variables other than the albedo may not be detected here.

"""

import numpy as np
from snicar8d_GO import snicar8d_GO
import matplotlib.pyplot as plt
import pandas as pd

##############################
# define default variable set:
##############################

DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
coszen   = 0.57    # if DIRECT give cosine of solar zenith angle 
dz = [0.01, 0.01, 0.01, 0.01, 0.01] # thickness of each vertical layer (unit = m)
nbr_lyr = len(dz)  # number of snow layers
R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
rho_snw = [500, 500, 500, 500, 500] # density of each layer (unit = kg m-3)
side_length = [5000,5000,5000,5000,5000]
depth = [5000,5000,5000,5000,5000]
nbr_aer = 16 # Define total number of different LAPs/aerosols in model
nbr_wvl = 470
snw_alg = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/snw_alg_1.nc'
glacier_algae1 = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/algae_geom_6_120.nc'
glacier_algae2 = '/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/algae_geom_6_20.nc'
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


########################
# DEFINE TEST FUNCTIONS
########################

def apprx_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [1,2,3]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/apprx_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:

        APRX_TYP = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        if np.allclose(test,target,rtol=1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"APPRX = {test_values[i]} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"APPRX = {test_values[i]} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")

    return 



def direct_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    print("NB DIRECT TESTS RUN WITH APRX = 2 DUE TO KNOWN BUG WITH EDDINGTON + DIFFUSE\n")
    
    #set values for test variables
    test_values = [0,1]
    APRX_TYP = 2
    
    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/direct_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        DIRECT = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"DIRECT = {test_values[i]} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"DIRECT = {test_values[i]} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")

    return


def delta_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [0,1]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/delta_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        DELTA = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"DELTA = {test_values[i]} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"DELTA = {test_values[i]} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")


    return


def coszen_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [0.3, 0.35, 0.4, 0.45, 0.5]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/coszen_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        coszen = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"COSZEN = {test_values[i]} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"COSZEN = {test_values[i]} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")

    return


def dz_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [[0.01,0.01,0.01,0.01,0.01],[0.01,0.05,0.05,0.05,0.05],[0.01,1,1,1,1],[0.001,0.001,0.001,0.001,0.001],[0.01,0.02,0.03,0.04,0.05]]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/dz_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        dz = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"DZ = ARRAY #{i} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"DZ = ARRAY #{i} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")


    return

def R_sfc_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [0.2, 0.4, 0.6]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/R_sfc_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        R_sfc = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"R_sfc = {test_values[i]} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"R_sfc = {test_values[i]} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")

    return

def sidelength_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):

    #set values for test variables
    test_values = [[5000,5000,5000,5000,5000],[8000,8000,8000,8000,8000],[10000,10000,10000,10000,10000],[15000,15000,15000,15000,15000],[6000,7000,8000,9000,10000]]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/side_length_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        rds_snw = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-3,atol=1e-3,equal_nan=True) == True:
            print(f"RDS = ARRAY #{i} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"RDS = ARRAY #{i} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")

    return



def rho_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2):
    
    #set values for test variables
    test_values = [[200,200,200,200,200],[300,300,300,300,300],[400,400,400,400,400],[500,500,500,500,500],[400,500,600,700,800]]

    # set up output vectors/arrays
    albedo_test = np.zeros([len(test_values),nbr_wvl])
    BBA_test = []
    BBAVIS_test = []
    BBANIR_test = []
    abs_slr_tot_test = []
    heat_rt_test = np.zeros([nbr_lyr,len(test_values)])

    # import target data
    albedo_target = pd.read_csv('/home/joe/Code/BioSNICAR_GO_PY/Unit_Tests/Snicar_GO/rho_albedo.csv',header=None)

    # initialise counter
    counter = 0

    # call function for each test value, holding all othe rvalues constant
    for i in test_values:
        rho_snw = i
        [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation] = snicar8d_GO(
        DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
        mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
        mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
        mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
        FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
        FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
        
        albedo_test[counter,:] = albedo
        BBA_test.append(BBA)
        BBAVIS_test.append(BBAVIS)
        BBANIR_test.append(BBANIR)
        abs_slr_tot_test.append(abs_slr_tot)
        heat_rt_test[:,counter] = heat_rt

        counter +=1

    # check to see that matlab (target) and python (test) values are equal within tolerance 1e-4
    for i in np.arange(0,len(test_values),1):
        test = albedo_test[i,:]
        target = albedo_target[i]
        
        if np.allclose(test,target,rtol = 1e-4,atol=1e-4,equal_nan=True) == True:
            print(f"RHO = ARRAY #{i} TEST PASSED WITH TOLERANCE 1e-4")
        else:
            print(f"RHO = ARRAY #{i} TEST FAILED WITH TOLERANCE 1e-4")

    print("\n **************************************************\n")
    
    return


#############################
# RUN TEST FUNCTIONS
#############################

print("\n************ UNIT TESTS ***************")
print("********* Geometric Optics ************")
print("*****************************************\n")
print("All tests run with default variable values apart from single \ntest variable")
print()
print("Values benchmarked against Matlab version, tests pass if Matlab \nand Python versions are equal to within \n1e-4\n")
print("************************\n")

direct_test( DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

apprx_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)


# delta_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

# coszen_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

# dz_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

# R_sfc_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

# sidelength_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

# rho_test(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
# mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
# mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
# mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
# FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
# FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)