
"""
##########################
SNICAR_LUT_builder.py
Joseph Cook, October 2019
##########################

This script calls SNICAR with a wide range of input variable values to populate a LUT with spectral albedos
for rapid model inversion (designd to couple with SNICAR_inverter.py).

Two functions are defined in this script. The first deals only with clean ice (i.e. no LAPs). SNICAR is called
with a range of grain side-lengths, depths and layer densities. The resulting spectral albedo is appended to a
ndarray. The spectral albedo is stored in the 4th dimension, with the first three being occupied by depth,
side-length and density. Since snicar is run layer-wise, each element of each dimension is an array of length n,
where n = number of layers. By default n=5. All variables other than depth, side length and density are held
constant at hard-coded values defined in-function.

The second function is identical to the first except that LAP mass mixing ratio is not held constant. This returns
a 4+n dimensional array where n = number of LAPs. In the case where the LAPs consist of a single dust and single 
algae, the result is a 6 dimensional array with the spectral albedo stored in the 6th dimension and the first five
dimensions occupied by depth, side-length, density, mass mixing ratio of dust, mass mixing ratio of algae. 

The resulting arrays are saved to disk as .npy files if save_LUT == True.

NB: The dirtyLUT builder talkes a long time to populate. For 8 depths, 8 side lengths, 8 densities, 9 dust1,
9 dust2 and 9 algae, 8*8*8*9*9 = 373248 runs are needed. Each run takes ~1.4 seconds (on i7-7700GHz processor),
meaning the entire LUT build takes 522547 seconds = 6.04 days. It would be helpful to find a way to parrelelize
this.

"""


import numpy as np
from snicar8d_GO import snicar8d_GO

##################################
# Set savepath and variable values
##################################

savepath = "/home/joe/Code/BioSNICAR_PY/Data/"

# set range of depths, side lengths, densities to iterate through: used in both clean & dirty functions
depths = [[3000,3000,3000,3000,3000],[5000,5000,5000,5000,5000],[7000,7000,7000,7000,7000],[9000,9000,9000,9000,9000],[12000,12000,12000,12000,12000],[15000,15000,15000,15000,15000]]
densities = [[400,400,400,400,400],[500,500,500,500,500],[600,600,600,600,600],[700,700,700,700,700],[800,800,800,800,800]]
wavelengths = np.arange(0.3,5,0.01)

# set dust and algae mass mixing ratios (ppb) to iterate through: only used in populate_array_dirty_ice()
algae = [[0,0,0,0,0],[5000,0,0,0,0],[10000,0,0,0,0],[25000,0,0,0,0],[50000,0,0,0,0],[75000,0,0,0,0],[100000,0,0,0,0],[125000,0,0,0,0],[150000,0,0,0,0],[200000,0,0,0,0]]
dust = [[0,0,0,0,0],[5000,0,0,0,0],[10000,0,0,0,0],[25000,0,0,0,0],[50000,0,0,0,0],[75000,0,0,0,0],[100000,0,0,0,0]]


##################
# DEFINE FUNCTIONS
##################


def populate_array_clean_ice(depths, side_lengths, densities, wavelengths, savepath, save_LUT = True):
    """
    param None
    return result: n-dimensional array with dim 1 = side length, dim 2 = depth, dim 3 = density, dim 4 = albedo

    This function populates a lookup table with snicar-predicted spectral albedos (0.3 - 5 micron) for a
    range of combinations of input variable values. This function considers only ice with no light-absorbing particles
    and only ice grains whose single scattering optical properties are generated using geometrical optics.

    """

    # set default variables (held constant - all variables other than depth, side-lengths and densities)
    DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
    coszen   = 0.57    # if DIRECT give cosine of solar zenith angle 
    dz = [0.001, 0.01, 0.02, 0.02, 0.02] # thickness of each vertical layer (unit = m)
    nbr_lyr = len(dz)  # number of snow layers
    R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
    nbr_aer = 16 # Define total number of different LAPs/aerosols in model
    nbr_wvl = 470
    snw_alg = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/snw_alg_1.nc'
    glacier_algae1 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_120.nc'
    glacier_algae2 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_20.nc'
    mss_cnc_soot1 = [0,0,0,0,0]    # uncoated black carbon
    mss_cnc_soot2 = [0,0,0,0,0]    # coated black carbon
    mss_cnc_dust1 = [0,0,0,0,0]    # global averside_lengths = [[3000,3000,3000,3000,3000],[5000,5000,5000,5000,5000],[7000,7000,7000,7000,7000],[9000,9000,9000,9000,9000],[12000,12000,12000,12000,12000],[15000,15000,15000,15000,15000]]age dust 1
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
    FILE_snw_alg  = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/snw_alg_1.nc"
    FILE_glacier_algae1 = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_120.nc" # Glacier algae
    FILE_glacier_algae2 = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_40.nc" # Glacier algae


    # set up result array
    result = np.zeros(shape = (len(depths),len(side_lengths),len(densities),len(wavelengths)))

    # loop through each combination of conditions and replace zeros with albedo at each wavelength
    # the result is a nd array with spectral albedo in the 4th dimension

    for i in np.arange(0,len(depths),1):
        for j in np.arange(0,len(side_lengths),1):
            for k in np.arange(0,len(densities),1):

                    depth = depths[i]
                    side_length = side_lengths[j]
                    rho_snw = densities[k]

                    wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot,heat_rt, total_insolation = snicar8d_GO(DIRECT, 
                    APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
                    mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
                    mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
                    mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
                    FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
                    FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
                    
                    result[i,j,k,:] = albedo
    
    
    if save_LUT:
        np.save(str(savepath+"clean_ice_LUT.npy"),result)

    return result



def populate_array_dirty_ice(depths, densities, wavelengths, dust, algae, savepath, save_LUT = True):

    # set default variables (held constant)
    DIRECT   = 1        # 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
    coszen   = 0.57    # if DIRECT give cosine of solar zenith angle 
    dz = [0.01, 0.01, 0.01, 0.01, 0.01] # thickness of each vertical layer (unit = m)
    nbr_lyr = len(dz)  # number of snow layers
    R_sfc = 0.15 # reflectance of undrlying surface - set across all wavelengths
    nbr_aer = 16 # Define total number of different LAPs/aerosols in model
    nbr_wvl = 470
    snw_alg = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/snw_alg_1.nc'
    glacier_algae1 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_120.nc'
    glacier_algae2 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_20.nc'
    mss_cnc_soot1 = [0,0,0,0,0]
    mss_cnc_soot2 = [0,0,0,0,0]    # coated black carbon
    mss_cnc_dust1 = [0,0,0,0,0]    # global average dust 1
    mss_cnc_dust2 = [0,0,0,0,0]    # global average dust 2
    mss_cnc_dust3 = [0,0,0,0,0]    # global average dust 3
    mss_cnc_dust4 = [0,0,0,0,0]    # global average dust 4
    mss_cnc_ash1 = [0,0,0,0,0]    # volcanic ash species 1
    mss_cnc_GRISdust2 = [0,0,0,0,0]    # GRIS dust 2 (Cook et al. 2019 HIGH)
    mss_cnc_GRISdust3 = [0,0,0,0,0]    # GRIS dust 3 (Cook et al. 2019 LOW)
    mss_cnc_GRISdustP1 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: low hematite)
    mss_cnc_GRISdustP2 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
    mss_cnc_GRISdustP3 = [0,0,0,0,0]  # GRIS dust 1 (Polashenki2015: median hematite)
    mss_cnc_snw_alg = [0,0,0,0,0]    # Snow Algae (spherical, C nivalis)
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
    FILE_snw_alg  = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/snw_alg_1.nc"
    FILE_glacier_algae1 = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_120.nc" # Glacier algae
    FILE_glacier_algae2 = "/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/algae_geom_6_40.nc" # Glacier algae

    # set up result array
    result = np.zeros(shape = (len(depths),len(densities),len(dust),len(algae),len(wavelengths)))

    # loop through each combination of conditions and replace zeros with albedo at each wavelength
    # the result is a nd array with spectral albedo in the 4th dimension

    for i in np.arange(0,len(depths),1):
            for j in np.arange(0,len(densities),1):
                for k in np.arange(0,len(dust),1):
                    for l in np.arange(0,len(algae),1):

                            depth = depths[i]
                            side_length = depths[i]
                            rho_snw = densities[j]
                            mss_cnc_GRISdust1 = dust[k]
                            mss_cnc_glacier_algae1 = algae[l]

                            wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot,heat_rt, total_insolation = snicar8d_GO(DIRECT, 
                            APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
                            mss_cnc_soot2, mss_cnc_dust1,mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
                            mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
                            mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
                            FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
                            FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)
                            
                            result[i,j,k,l,:] = albedo
    
    if save_LUT:
        np.save("/home/joe/Desktop/dirty_ice_LUT_BIG.npy",result)
    
    return result

#######################
#### RUN FUNCTIONS
#######################

#clean_result = populate_array_clean_ice(depths, side_lengths, densities, wavelengths, savepath, save_LUT = False)
dirty_result = populate_array_dirty_ice(depths, densities, wavelengths, dust, algae, savepath, save_LUT = True)
