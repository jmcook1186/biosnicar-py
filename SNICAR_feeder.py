def snicar_feeder(MIE, GO, dir_base, rf_ice, incoming_i, DIRECT, layer_type,\
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
    FILE_snw_alg, FILE_glacier_algae):


    """
    This script takes the user defined inputs from the driver script and calculates te relevant
    tau (optical thickness), SSA (single scattering albedo) and g (asymmetry parameter) for each 
    vertical layer to send to the radiative transfer solver.

    There are two options for the radiative transfer solver A) the fast tridiaginal matrix method 
    of Toon et al. (1989) which is very efficient but limited to granular layers
    or B) the adding-doubling method that can include bubbly ice layers and fresnel
    reflections anywhere in the column.

    For Mie calculations, this script makes the necessary ajustments for nonspherical grain shapes
    using the method of He et al. (2016).

    The script calls out to one of two radiative transfer solver scripts: adding_doubling_solver.py
    or two_stream_solver.py.

    """


    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    from IceOptical_Model.mie_coated_water_spheres import miecoated_driver
    from Toon_RT_solver import toon_solver
    from adding_doubling_solver import adding_doubling_solver

    # set working directory (location of netcdf library)
    dir_mie_files = "Data/Mie_files/"
    dir_GO_files = 'Data/GO_files/'

    # retrieve wavelength from arbitrary choice of netcdf file
    temp = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/ice_Wrn84/ice_Wrn84_0500.nc"))
    wvl = np.array(temp['wvl'].values)
    wvl = wvl*1e6
    nbr_wvl = len(wvl)

    # set reflectance of underlying surface
    R_sfc = [R_sfc for _ in range(nbr_wvl)]
    R_sfc = np.array(R_sfc)

    # Incoming Irradiance
    # calc cosine of solar zenith (radians)
    mu_not = np.cos(solzen * (np.pi / 180)) # convert radians if required
    
    print("\ncosine of solar zenith = ", mu_not)
    
    flx_slr = []
    
    if DIRECT:

        if incoming_i == 0:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_mlw_clr.nc"))
            print("atmospheric profile = mid-lat winter")
        elif incoming_i == 1:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_mls_clr.nc"))
            print("atmospheric profile = mid-lat summer")
        elif incoming_i == 2:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_saw_clr.nc"))
            print("atmospheric profile = sub-Arctic winter")
        elif incoming_i == 3:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_sas_clr.nc"))
            print("atmospheric profile = sub-Arctic summer")
        elif incoming_i == 4:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_smm_clr.nc"))
            print("atmospheric profile = Summit Station")
        elif incoming_i == 5:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_hmn_clr.nc"))
            print("atmospheric profile = High Mountain")  
        elif incoming_i == 6:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_toa_clr.nc"))
            print("atmospheric profile = top-of-atmosphere")

        else:
            raise ValueError ("Invalid choice of atmospheric profile")
        
        flx_slr = Incoming_file['flx_dwn_sfc'].values
        flx_slr[flx_slr<=0]=1e-30
        Fs = flx_slr / (mu_not * np.pi)
        Fd = np.zeros(nbr_wvl)

    else:

        if incoming_i == 0:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_mlw_cld.nc"))
        elif incoming_i == 1:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_mls_cld.nc"))
        elif incoming_i == 2:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_saw_cld.nc"))
        elif incoming_i == 3:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_sas_cld.nc"))
        elif incoming_i == 4:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_smm_cld.nc"))
        elif incoming_i == 5:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_hmn_cld.nc"))   
        elif incoming_i == 6:
            Incoming_file = xr.open_dataset(str(dir_base+dir_mie_files+"/480band/fsds/swnb_480bnd_toa_cld.nc"))

        else:
            raise ValueError ("Invalid choice of atmospheric profile")     

        flx_slr = Incoming_file['flx_dwn_sfc'].values
        
        flx_slr[flx_slr<=0]=1e-30

        Fd = [flx_slr[i]/mu_not*np.pi for i in range(nbr_wvl)]
        Fs = np.zeros(nbr_wvl)


    # Read in ice optical properties
    # set string stubs for reading in ice optical files

    if rf_ice == 0:
        fl_stb1 = '480band/ice_Wrn84/ice_Wrn84_'
        print("Using Warren 84 refractive index")
    elif rf_ice == 1:
        fl_stb1 = '480band/ice_Wrn08/ice_Wrn08_'
        print("Using Warren 08 refractive index")
    elif rf_ice == 2:
        fl_stb1 = '480band/ice_Pic16/ice_Pic16_'
        print("Using Picard 16 refractive index")

    fl_stb2 = ".nc"


    # set up empty arrays
    SSA_snw = np.empty([nbr_lyr, nbr_wvl])
    MAC_snw = np.empty([nbr_lyr, nbr_wvl])
    g_snw = np.empty([nbr_lyr, nbr_wvl])
    abs_cff_mss_ice = np.empty([nbr_wvl])

    for i in np.arange(0,nbr_lyr,1):

        if layer_type[i] == 0: # (granular layer)

            if rds_snw[i] ==0:

                raise ValueError("ERROR: ICE GRAIN RADIUS SET TO ZERO")

            else:

                if MIE:
                    s1 = "{}".format(str(rds_snw[i]))
                    s1 = s1.rjust(4,'0')
                    FILE_ice = str(dir_base + dir_mie_files + fl_stb1 + s1 + fl_stb2)
                    print("\nLayer: {}".format(i))
                    print("Using Mie mode: spheres with radius = {}".format(s1))

                elif GO:

                    s1 = "{}".format(str(side_length[i]))
                    s1 = s1.rjust(4,'0')
                    s2 = "{}".format(str(depth[i]))
                    FILE_ice = str(dir_base + dir_GO_files + fl_stb1 + s1 + "_" + s2 + fl_stb2)
                    print("\nLayer: {}".format(i))
                    print("Using geometric optics mode: hex column with side length = {}, length = {}".format(s1,s2))

            # read in single scattering albedo, MAC and g for ice crystals in each layer,
            # optional with coated liquid water spheres
            if rwater[i] > rds_snw[i]:

                # water coating code currently disabled
                raise ValueError("SORRY, water coatings are not yet functional until the code has been updated to deal with wl's down to 200 nm")               
                
                #    fn_ice = dir_base + "/Data/rfidx_ice.nc"

                #    fn_water = dir_base + "Data/Refractive_Index_Liquid_Water_Segelstein_1981.csv"
                #     res = miecoated_driver(rice=rds_snw[i], rwater=rwater[i], fn_ice=fn_ice, fn_water=fn_water, wvl=wvl)
                #     SSA_snw[i, :] = res["ssa"]
                #     g_snw[i, :] = res["asymmetry"]

                #     with xr.open_dataset(FILE_ice) as temp:
                #         ext_cff_mss = temp['ext_cff_mss'].values
                #         MAC_snw[i, :] = ext_cff_mss

            else:

                with xr.open_dataset(FILE_ice) as temp:
                    
                    SSA = temp['ss_alb'].values
                    SSA_snw[i,:] = SSA

                    ext_cff_mss = temp['ext_cff_mss'].values
                    MAC_snw[i,:] = ext_cff_mss

                    asm_prm = temp['asm_prm'].values
                    g_snw[i,:] = asm_prm
        

            ###############################################################
            ## NEW GRAIN SHAPE ROUTINE
            ################################################################

            # Constants for aspherical ice particles
            # g_snw asymmetry factor parameterization coefficients (6 bands) from
            # Table 3 & Eqs. 6-7 in He et al. (2017)
            # assume same values for 4-5 um band, which leads to very small biases (<3%)
            
            g_wvl = np.array([0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00]) # wavelength (um) division point
            g_wvl_center = np.array(g_wvl[1:8])/2 + np.array(g_wvl[0:7])/2  # center point for wavelength band
            g_b0 = np.array([9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01])
            g_b1 = np.array([5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02])
            g_b2 = np.array([-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04])
            
            # Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007
            g_F07_c2 = np.array([1.349959e-1,1.115697e-1,9.853958e-2,5.557793e-2,-1.233493e-1,0.0,0.0])
            g_F07_c1 = np.array([-3.987320e-1,-3.723287e-1,-3.924784e-1,-3.259404e-1,4.429054e-2,-1.726586e-1,-1.726586e-1])
            g_F07_c0 = np.array([7.938904e-1,8.030084e-1,8.513932e-1,8.692241e-1,7.085850e-1,6.412701e-1,6.412701e-1])
            g_F07_p2 = np.array([3.165543e-3,2.014810e-3,1.780838e-3,6.987734e-4,-1.882932e-2,-2.277872e-2,-2.277872e-2])
            g_F07_p1 = np.array([1.140557e-1,1.143152e-1,1.143814e-1,1.071238e-1,1.353873e-1,1.914431e-1,1.914431e-1])
            g_F07_p0 = np.array([5.292852e-1,5.425909e-1,5.601598e-1,6.023407e-1,6.473899e-1,4.634944e-1,4.634944e-1])
            fs_hex = 0.788 # shape factor for hexagonal plate (reference)
            

            for i in np.arange(0,nbr_lyr,1):
                
                if layer_type[i] == 0:
                    
                    if snw_shp[i] == 0: # snow

                        pass # if layer type is for spheres, no changes required

                    
                    elif snw_shp[i] == 1: # 1 = spheroid, He et al. (2017) parameterization
                        
                        diam_ice = 2.0 * rds_snw[i] # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:

                            fs_sphd = 0.929 # default shape factor for spheroid; He et al. (2017), Table 1
                        
                        else:
                            
                            fs_sphd = shp_fctr[i] # if shp_factor not 0, then use user-defined value
                        
                        
                        if snw_ar[i] == 0:

                            AR_tmp = 0.5 # default aspect ratio for spheroid; He et al. (2017), Table 1
                        
                        else:

                            AR_tmp = snw_ar[i]
                        
                        g_snw_Cg_tmp = g_b0 * (fs_sphd/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * AR_tmp**2 # Eqn. 3.1 in Fu (2007)


                    elif snw_shp[i] == 2: # 3=hexagonal plate, He et al. 2017 parameterization

                        diam_ice = 2.0 * rds_snw[i] # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:
                            
                            fs_hex0 = 0.788 # default shape factor for hexagonal plates; He et al. (2017), Table 1
                        
                        else:

                            fs_hex0 = shp_fctr[i]
                        

                        if snw_ar[i] == 0:

                            AR_tmp = 2.5 # default aspect ratio for hexagonal plate; He et al. (2017), Table 1
                        
                        else:
                        
                            AR_tmp = snw_ar[i]
                                
                        g_snw_Cg_tmp = g_b0 * (fs_hex0/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2   # Eqn. 3.3 in Fu (2007)
                        

                    elif snw_shp[i] == 3: # 4=koch snowflake, He et al. (2017) parameterization

                        diam_ice = 2.0 * rds_snw[i] / 0.544 # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:
                        
                            fs_koch = 0.712 # default shape factor for koch snowflake; He et al. (2017), Table 1
                        
                        else:
                            
                            fs_koch = shp_fctr[i]
                        

                        if snw_ar[i] == 0:

                            AR_tmp = 2.5 # default aspect ratio for koch snowflake; He et al. (2017), Table 1
                        
                        else:

                            AR_tmp = snw_ar[i]
                        
                        g_snw_Cg_tmp = g_b0 * (fs_koch/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2  # Eqn. 3.3 in Fu (2007)
                        
                    
                    if snw_shp[i] > 0:

                        from scipy.interpolate import pchip
                        # 6 wavelength bands for g_snw to be interpolated into 480-bands of SNICAR
                        # shape-preserving piecewise interpolation into 480-bands
                        g_Cg_intp = pchip(g_wvl_center,g_snw_Cg_tmp)(wvl)
                        gg_F07_intp = pchip(g_wvl_center,gg_snw_F07_tmp)(wvl)
                        g_snw_F07 = gg_F07_intp + (1.0 - gg_F07_intp) / SSA_snw[i,:] / 2 # Eq.2.2 in Fu (2007)
                        g_snw[i,:] = g_snw_F07 * g_Cg_intp # Eq.6, He et al. (2017)
                        g_snw[i,371:470] = g_snw[i,370] # assume same values for 4-5 um band, with very small biases (<3%)
                    
                    g_snw[g_snw < 0] = 0.01
                    g_snw[g_snw > 0.99] = 0.99 # avoid unreasonable values (so far only occur in large-size spheroid cases)


        else: # else correspondng to if layer_type == 1
            
            rd = "{}".format(rds_snw[i])
            rd = rd.rjust(4,"0")
            refidx_file = xr.open_dataset('/home/joe/Code/BioSNICAR_GO_PY/Data/rfidx_ice.nc')
           
            if rf_ice == 0:
                refidx_re = refidx_file['re_Wrn84'].values
                refidx_im = refidx_file['im_Wrn84'].values 

            elif rf_ice == 1:
                refidx_re = refidx_file['re_Wrn08'].values
                refidx_im = refidx_file['im_Wrn08'].values 

            elif rf_ice == 2:
                refidx_re = refidx_file['re_Pic16'].values
                refidx_im = refidx_file['im_Pic16'].values 

            FILE_ice = '/home/joe/Code/BioSNICAR_GO_PY/Data/bubbly_ice_files/bbl_{}.nc'.format(rd)
            file = xr.open_dataset(FILE_ice)
            sca_cff_vlm = file['sca_cff_vlm'].values # scattering cross section unit per volume of bubble
            g_snw[i,:] = file['asm_prm'].values
            abs_cff_mss_ice[:] = ((4 * np.pi * refidx_im) / (wvl * 1e-6))/917
            vlm_frac_air = (917 - rho_snw[i]) / 917
            MAC_snw[i,:] = ((sca_cff_vlm * vlm_frac_air) /917) + abs_cff_mss_ice
            SSA_snw[i,:] = ((sca_cff_vlm * vlm_frac_air) /917) / MAC_snw[i,:]


    
    ###############################################################
    ##############################################################

    # open netcdf files
    FILE_soot1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_soot1))
    FILE_soot2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_soot2))
    FILE_brwnC1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_brwnC1))
    FILE_brwnC2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_brwnC2))
    FILE_dust1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_dust1))
    FILE_dust2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_dust2))
    FILE_dust3 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_dust3))
    FILE_dust4 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_dust4))
    FILE_ash1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_ash1))
    FILE_ash2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_ash2))
    FILE_ash3 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_ash3))
    FILE_ash4 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_ash4))
    FILE_ash5 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_ash5))
    FILE_Skiles_dust1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_Skiles_dust1))
    FILE_Skiles_dust2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_Skiles_dust2))
    FILE_Skiles_dust3 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_Skiles_dust3))
    FILE_Skiles_dust4 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_Skiles_dust4))
    FILE_Skiles_dust5 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_Skiles_dust5))
    FILE_GreenlandCentral1 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_GreenlandCentral1))
    FILE_GreenlandCentral2 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_GreenlandCentral2))
    FILE_GreenlandCentral3 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_GreenlandCentral3))
    FILE_GreenlandCentral4 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_GreenlandCentral4))
    FILE_GreenlandCentral5 = xr.open_dataset(str(dir_base + dir_mie_files+ FILE_GreenlandCentral5))
    FILE_snw_alg = xr.open_dataset(str(dir_base + dir_mie_files+FILE_snw_alg))
    FILE_glacier_algae = xr.open_dataset(str(dir_base + dir_mie_files+FILE_glacier_algae))

    # read in aerosol optical properties
    SSAaer = np.zeros([nbr_aer,nbr_wvl])
    SSAaer[0,:] = FILE_soot1['ss_alb'].values
    SSAaer[1,:] = FILE_soot2['ss_alb'].values
    SSAaer[2,:] = FILE_brwnC1['ss_alb'].values
    SSAaer[3,:] = FILE_brwnC2['ss_alb'].values
    SSAaer[4,:] = FILE_dust1['ss_alb'].values
    SSAaer[5,:] = FILE_dust2['ss_alb'].values
    SSAaer[6,:] = FILE_dust3['ss_alb'].values
    SSAaer[7,:] = FILE_dust4['ss_alb'].values
    SSAaer[8,:] = FILE_ash1['ss_alb'].values
    SSAaer[9,:] = FILE_ash2['ss_alb'].values
    SSAaer[10,:] = FILE_ash3['ss_alb'].values
    SSAaer[11,:] = FILE_ash4['ss_alb'].values
    SSAaer[12,:] = FILE_ash5['ss_alb'].values
    SSAaer[13,:] = FILE_Skiles_dust1['ss_alb'].values
    SSAaer[14,:] = FILE_Skiles_dust2['ss_alb'].values
    SSAaer[15,:] = FILE_Skiles_dust3['ss_alb'].values
    SSAaer[16,:] = FILE_Skiles_dust4['ss_alb'].values
    SSAaer[17,:] = FILE_Skiles_dust5['ss_alb'].values
    SSAaer[18,:] = FILE_GreenlandCentral1['ss_alb'].values
    SSAaer[19,:] = FILE_GreenlandCentral2['ss_alb'].values
    SSAaer[20,:] = FILE_GreenlandCentral3['ss_alb'].values
    SSAaer[21,:] = FILE_GreenlandCentral4['ss_alb'].values
    SSAaer[22,:] = FILE_GreenlandCentral5['ss_alb'].values
    SSAaer[23,:] = FILE_snw_alg['ss_alb'].values
    SSAaer[24,:] = FILE_glacier_algae['ss_alb'].values

    MACaer = np.zeros([nbr_aer, nbr_wvl])

    MACaer[0,:] = FILE_soot1['ext_cff_mss'].values 
    MACaer[1,:] = FILE_soot2['ext_cff_mss_ncl'].values # normalised to mass of carbon, not carbon plus coating
    MACaer[2,:] = FILE_brwnC1['ext_cff_mss'].values
    MACaer[3,:] = FILE_brwnC2['ext_cff_mss_ncl'].values
    MACaer[4,:] = FILE_dust1['ext_cff_mss'].values
    MACaer[5,:] = FILE_dust2['ext_cff_mss'].values
    MACaer[6,:] = FILE_dust3['ext_cff_mss'].values
    MACaer[7,:] = FILE_dust4['ext_cff_mss'].values
    MACaer[8,:] = FILE_ash1['ext_cff_mss'].values
    MACaer[9,:] = FILE_ash2['ext_cff_mss'].values
    MACaer[10,:] = FILE_ash3['ext_cff_mss'].values
    MACaer[11,:] = FILE_ash4['ext_cff_mss'].values
    MACaer[12,:] = FILE_ash5['ext_cff_mss'].values
    MACaer[13,:] = FILE_Skiles_dust1['ext_cff_mss'].values
    MACaer[14,:] = FILE_Skiles_dust2['ext_cff_mss'].values
    MACaer[15,:] = FILE_Skiles_dust3['ext_cff_mss'].values
    MACaer[16,:] = FILE_Skiles_dust4['ext_cff_mss'].values
    MACaer[17,:] = FILE_Skiles_dust5['ext_cff_mss'].values
    MACaer[18,:] = FILE_GreenlandCentral1['ext_cff_mss'].values
    MACaer[19,:] = FILE_GreenlandCentral2['ext_cff_mss'].values
    MACaer[20,:] = FILE_GreenlandCentral3['ext_cff_mss'].values
    MACaer[21,:] = FILE_GreenlandCentral4['ext_cff_mss'].values
    MACaer[22,:] = FILE_GreenlandCentral5['ext_cff_mss'].values
    MACaer[23,:] = FILE_snw_alg['ext_cff_mss'].values
    MACaer[24,:] = FILE_glacier_algae['ext_cff_mss'].values

    Gaer = np.zeros([nbr_aer,nbr_wvl])

    Gaer[0,:] = FILE_soot1['asm_prm'].values
    Gaer[1,:] = FILE_soot2['asm_prm'].values
    Gaer[2,:] = FILE_brwnC1['asm_prm'].values
    Gaer[3,:] = FILE_brwnC2['asm_prm'].values
    Gaer[4,:] = FILE_dust1['asm_prm'].values
    Gaer[5,:] = FILE_dust2['asm_prm'].values
    Gaer[6,:] = FILE_dust3['asm_prm'].values
    Gaer[7,:] = FILE_dust4['asm_prm'].values
    Gaer[8,:] = FILE_ash1['asm_prm'].values
    Gaer[9,:] = FILE_ash2['asm_prm'].values
    Gaer[10,:] = FILE_ash3['asm_prm'].values
    Gaer[11,:] = FILE_ash4['asm_prm'].values
    Gaer[12,:] = FILE_ash5['asm_prm'].values
    Gaer[13,:] = FILE_Skiles_dust1['asm_prm'].values
    Gaer[14,:] = FILE_Skiles_dust2['asm_prm'].values
    Gaer[15,:] = FILE_Skiles_dust3['asm_prm'].values
    Gaer[16,:] = FILE_Skiles_dust4['asm_prm'].values
    Gaer[17,:] = FILE_Skiles_dust5['asm_prm'].values
    Gaer[18,:] = FILE_GreenlandCentral1['asm_prm'].values
    Gaer[19,:] = FILE_GreenlandCentral2['asm_prm'].values
    Gaer[20,:] = FILE_GreenlandCentral3['asm_prm'].values
    Gaer[21,:] = FILE_GreenlandCentral4['asm_prm'].values
    Gaer[22,:] = FILE_GreenlandCentral5['asm_prm'].values
    Gaer[23,:] = FILE_snw_alg['asm_prm'].values
    Gaer[24,:] = FILE_glacier_algae['asm_prm'].values


    # load mass concentrations per layer into numpy array (one row per layer, one column per umpurity)
    # and convert to kg/kg unit

    MSSaer = np.zeros([nbr_lyr, nbr_aer])
    MSSaer[0:nbr_lyr,0] = mss_cnc_soot1
    MSSaer[0:nbr_lyr,1] = mss_cnc_soot2
    MSSaer[0:nbr_lyr,2] = mss_cnc_brwnC1
    MSSaer[0:nbr_lyr,3] = mss_cnc_brwnC2
    MSSaer[0:nbr_lyr,4] = mss_cnc_dust1
    MSSaer[0:nbr_lyr,5] = mss_cnc_dust2
    MSSaer[0:nbr_lyr,6] = mss_cnc_dust3
    MSSaer[0:nbr_lyr,7] = mss_cnc_dust4
    MSSaer[0:nbr_lyr,8] = mss_cnc_ash1
    MSSaer[0:nbr_lyr,9] = mss_cnc_ash2
    MSSaer[0:nbr_lyr,10] = mss_cnc_ash3
    MSSaer[0:nbr_lyr,11] = mss_cnc_ash4
    MSSaer[0:nbr_lyr,12] = mss_cnc_ash5
    MSSaer[0:nbr_lyr,13] = mss_cnc_Skiles_dust1
    MSSaer[0:nbr_lyr,14] = mss_cnc_Skiles_dust2
    MSSaer[0:nbr_lyr,15] = mss_cnc_Skiles_dust3
    MSSaer[0:nbr_lyr,16] = mss_cnc_Skiles_dust4
    MSSaer[0:nbr_lyr,17] = mss_cnc_Skiles_dust5
    MSSaer[0:nbr_lyr,18] = mss_cnc_GreenlandCentral1
    MSSaer[0:nbr_lyr,19] = mss_cnc_GreenlandCentral2
    MSSaer[0:nbr_lyr,20] = mss_cnc_GreenlandCentral3
    MSSaer[0:nbr_lyr,21] = mss_cnc_GreenlandCentral4
    MSSaer[0:nbr_lyr,22] = mss_cnc_GreenlandCentral5
    MSSaer[0:nbr_lyr,23] = mss_cnc_snw_alg
    MSSaer[0:nbr_lyr,24] = mss_cnc_glacier_algae

    MSSaer = MSSaer*1e-9


    #####################################
    # Begin solving Radiative Transfer
    #####################################

    """
    #1. Calculate effective tau (optical depth), SSA (single scattering albedo) and 
    # g (assymetry parameter) for the ice + impurities mixture.

    # SSA and g for the individual components has been calculated using Mie theory and
    # stored in a netcdf file. Here, these values are combined to give an overall
    # SSA and g for the ice + impurity mixture
    
    """

    # initialize arrays
    g_sum = np.zeros([nbr_lyr, nbr_wvl])
    SSA_sum = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau = np.zeros([nbr_lyr, nbr_wvl])
    SSA = np.zeros([nbr_lyr, nbr_wvl])
    g = np.zeros([nbr_lyr, nbr_wvl])
    L_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_aer = np.zeros([nbr_lyr, nbr_aer, nbr_wvl])
    tau_sum = np.zeros([nbr_lyr, nbr_wvl])
    SSA_sum = np.zeros([nbr_lyr, nbr_wvl])
    L_snw = np.zeros(nbr_lyr)
    tau_snw = np.zeros([nbr_lyr,nbr_wvl])


    # for each layer, the layer mass (L) is density * layer thickness
    # for each layer the optical depth is the layer mass * the mass extinction coefficient
    # first for the ice in each layer
    
    for i in range(nbr_lyr):

        L_snw[i] = rho_snw[i] * dz[i]
        tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]

    # then for the LAPs in each layer
    for i in range(nbr_lyr):
        for j in range(nbr_aer):

            L_aer[i, j, :] = L_snw[i] * MSSaer[i, j]
            tau_aer[i, j, :] = L_aer[i, j, :] * MACaer[j, :]

            tau_sum = tau_sum + tau_aer[i, j, :]
            SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
            g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])


    # finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
    for i in range(nbr_lyr):
        
        tau[i,:] = tau_sum[i,:] + tau_snw[i,:]
        SSA[i,:] = (1 / tau[i,:]) * (SSA_sum[i,:] + SSA_snw[i,:] * tau_snw[i,:])
        g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i,:] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))

    # just in case any unrealistic values arise (none detected so far)
    SSA[SSA<=0]=0.00000001
    SSA[SSA>=1]=0.99999999
    g[g<=0]=0.00001
    g[g>=1]=0.99999

    # CALL RT SOLVER (TOON  = TOON ET AL, TRIDIAGONAL MATRIX METHOD; 
    # ADD_DOUBLE = ADDING-DOUBLING METHOD)
   
    if TOON: 

        wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_vis_tot, heat_rt = \
            toon_solver(APRX_TYP, DELTA, tau, g, SSA, mu_not, nbr_lyr, nbr_wvl, R_sfc, wvl, Fs, Fd,\
            L_snw, flx_slr)


    if ADD_DOUBLE:

        wvl, flx_dwn_spc, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt = \
            adding_doubling_solver(rf_ice, APRX_TYP, DELTA, layer_type, tau, g, SSA, mu_not, nbr_lyr, nbr_wvl, R_sfc, wvl, Fs, Fd,\
            L_snw, flx_slr, DIRECT, rds_snw)


    return wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt