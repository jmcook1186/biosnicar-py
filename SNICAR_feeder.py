def snicar_feeder(inputs):
    

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
    import random
    import os
    import collections as c
    
    # load variables from input table
    dir_base=inputs.dir_base
    TOON=inputs.TOON
    ADD_DOUBLE=inputs.ADD_DOUBLE
    nbr_lyr = inputs.nbr_lyr
    nbr_aer = inputs.nbr_aer
    layer_type = inputs.layer_type
    DIRECT=inputs.DIRECT
    incoming_i=inputs.incoming_i
    solzen=inputs.solzen
    grain_rds=inputs.grain_rds
    grain_shp=inputs.grain_shp
    shp_fctr=inputs.shp_fctr
    grain_ar=inputs.grain_ar
    side_length=inputs.side_length
    depth=inputs.depth
    rf_ice=inputs.rf_ice
    rwater=inputs.rwater
    rho_layers=inputs.rho_layers
    dz=inputs.dz
    FILE_brwnC2=inputs.FILE_brwnC2
    FILE_soot2=inputs.FILE_soot2
    GA_units = inputs.GA_units
    Cfactor = inputs.Cfactor
    
    
    files = [inputs.FILE_soot1,\
    inputs.FILE_soot2, inputs.FILE_brwnC1, inputs.FILE_brwnC2, inputs.FILE_dust1, inputs.FILE_dust2, inputs.FILE_dust3, inputs.FILE_dust4, inputs.FILE_dust5,\
    inputs.FILE_ash1, inputs.FILE_ash2, inputs.FILE_ash3, inputs.FILE_ash4, inputs.FILE_ash5, inputs.FILE_ash_st_helens, inputs.FILE_Skiles_dust1, inputs.FILE_Skiles_dust2,\
    inputs.FILE_Skiles_dust3, inputs.FILE_Skiles_dust4, inputs.FILE_Skiles_dust5, inputs.FILE_GreenlandCentral1,\
    inputs.FILE_GreenlandCentral2, inputs.FILE_GreenlandCentral3, inputs.FILE_GreenlandCentral4, inputs.FILE_GreenlandCentral5,\
    inputs.FILE_Cook_Greenland_dust_L, inputs.FILE_Cook_Greenland_dust_C, inputs.FILE_Cook_Greenland_dust_H,\
    inputs.FILE_snw_alg, inputs.FILE_glacier_algae]
        
    mass_concentrations = [inputs.mss_cnc_soot1, inputs.mss_cnc_soot2, inputs.mss_cnc_brwnC1, inputs.mss_cnc_brwnC2, inputs.mss_cnc_dust1,\
    inputs.mss_cnc_dust2, inputs.mss_cnc_dust3, inputs.mss_cnc_dust4, inputs.mss_cnc_dust5, inputs.mss_cnc_ash1, inputs.mss_cnc_ash2,\
    inputs.mss_cnc_ash3, inputs.mss_cnc_ash4, inputs.mss_cnc_ash5, inputs.mss_cnc_ash_st_helens, inputs.mss_cnc_Skiles_dust1, inputs.mss_cnc_Skiles_dust2,\
    inputs.mss_cnc_Skiles_dust3, inputs.mss_cnc_Skiles_dust4, inputs.mss_cnc_Skiles_dust5, inputs.mss_cnc_GreenlandCentral1,\
    inputs.mss_cnc_GreenlandCentral2, inputs.mss_cnc_GreenlandCentral3, inputs.mss_cnc_GreenlandCentral4,\
    inputs.mss_cnc_GreenlandCentral5, inputs.mss_cnc_Cook_Greenland_dust_L, inputs.mss_cnc_Cook_Greenland_dust_C,\
    inputs.mss_cnc_Cook_Greenland_dust_H, inputs.mss_cnc_snw_alg, inputs.mss_cnc_glacier_algae]
        
    # working directories 
    dir_mie_ice_files = str(dir_base + 'Data/Mie_files/480band/') # directory with folders ice_Pic16, ice_Wrn08 and ice_Wrn84 with optical properties calculated with Mie theory
    dir_go_ice_files = str(dir_base + 'Data/GO_files/480band/') # idem for ice OPs calculated with Geometric optics
    dir_mie_lap_files = str(dir_base + 'Data/Mie_files/480band/lap/') # directory with folders ice_Pic16, ice_Wrn08 and ice_Wrn84 with optical properties calculated with Mie theory
    dir_bubbly_ice = str(dir_base + 'Data/bubbly_ice_files/')
    dir_fsds = str(dir_base + 'Data/Mie_files/480band/fsds/')
    dir_RI_ice = str(dir_base + 'Data/') 

    # retrieve nbr wvl, aer, layers and layer types 
    temp = xr.open_dataset(str(dir_mie_lap_files+random.choice(os.listdir(dir_mie_lap_files))))
    wvl = np.array(temp['wvl'].values)
    wvl = wvl*1e6
    nbr_wvl = len(wvl)
    inputs.nbr_wvl = nbr_wvl
    inputs.wvl = wvl

    # load incoming irradiance
    # calc cosine of solar zenith (radians)
    mu_not = np.cos(solzen * (np.pi / 180)) # convert radians if required
    inputs.mu_not = mu_not
    
    print("\ncosine of solar zenith = ", mu_not)
    
    flx_slr = []
    
    if DIRECT:
        
        
        coszen = str('SZA'+str(solzen).rjust(2,'0'))

        if incoming_i == 0:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mlw_clr_"+coszen+".nc")) 
            print("atmospheric profile = mid-lat winter")
        elif incoming_i == 1:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mls_clr_"+coszen+".nc"))
            print("atmospheric profile = mid-lat summer")
        elif incoming_i == 2:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_saw_clr_"+coszen+".nc"))
            print("atmospheric profile = sub-Arctic winter")
        elif incoming_i == 3:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_sas_clr_"+coszen+".nc"))
            print("atmospheric profile = sub-Arctic summer")
        elif incoming_i == 4:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_smm_clr_"+coszen+".nc"))
            print("atmospheric profile = Summit Station")
        elif incoming_i == 5:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_hmn_clr_"+coszen+".nc"))
            print("atmospheric profile = High Mountain")  
        elif incoming_i == 6:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_toa_clr.nc"))
            print("atmospheric profile = top-of-atmosphere")

        else:
            raise ValueError ("Invalid choice of atmospheric profile")
        
        flx_slr = Incoming_file['flx_dwn_sfc'].values #flx_dwn_sfc is the spectral irradiance in W m-2 and is pre-calculated (flx_frc_sfc*flx_bb_sfc in original code)
        flx_slr[flx_slr<=0]=1e-30
        inputs.flx_slr=flx_slr
        inputs.Fs = flx_slr / (mu_not * np.pi)
        inputs.Fd = np.zeros(nbr_wvl)

    else:

        if incoming_i == 0:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mlw_cld.nc"))
        elif incoming_i == 1:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_mls_cld.nc"))
        elif incoming_i == 2:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_saw_cld.nc"))
        elif incoming_i == 3:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_sas_cld.nc"))
        elif incoming_i == 4:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_smm_cld.nc"))
        elif incoming_i == 5:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_hmn_cld.nc"))   
        elif incoming_i == 6:
            Incoming_file = xr.open_dataset(str(dir_fsds + "swnb_480bnd_toa_cld.nc"))

        else:
            raise ValueError ("Invalid choice of atmospheric profile")     

        flx_slr = Incoming_file['flx_dwn_sfc'].values
        
        flx_slr[flx_slr<=0]=1e-30

        inputs.Fd = [flx_slr[i]/mu_not*np.pi for i in range(nbr_wvl)]
        inputs.Fs = np.zeros(nbr_wvl)


    ###################################################
    # Read in ice optical properties
    ###################################################
    # set up empty arrays
    SSA_snw = np.empty([nbr_lyr, nbr_wvl])
    MAC_snw = np.empty([nbr_lyr, nbr_wvl])
    g_snw = np.empty([nbr_lyr, nbr_wvl])
    abs_cff_mss_ice = np.empty([nbr_wvl])

    for i in np.arange(0,nbr_lyr,1):

        if layer_type[i] == 0: # (granular layer)

            if grain_rds[i] == 0:

                raise ValueError("ERROR: ICE GRAIN RADIUS SET TO ZERO")

            else:
                
                if grain_shp[i] == 4: # if large hexaginal prisms (geometric optics calcs)
                    print("Using hexagonal column with side length = {}, length = {}".format(str(side_length[i]).rjust(4,'0'),str(depth[i])))
                  
                    if rf_ice ==0:
                        dir_OP = str(dir_go_ice_files+'ice_Wrn84/ice_Wrn84_')
                        print("Using Warren 84 refractive index")
                    elif rf_ice == 1:
                        dir_OP = str(dir_go_ice_files+'ice_Wrn08/ice_Wrn08_')
                        print("Using Warren 08 refractive index")
                    elif rf_ice == 2:
                        dir_OP = str(dir_go_ice_files+'ice_Pic16/ice_Pic16_')
                        print("Using Picard 16 refractive index")           

                    FILE_ice = str(dir_OP + '{}_{}.nc'.format(str(side_length[i]).rjust(4,'0'), str(depth[i])))
                    print("\nLayer: {}".format(i))
    


                elif grain_shp[i] < 4:

                    if rf_ice == 0:
                        dir_OP = 'ice_Wrn84/ice_Wrn84'
                        print("Using Warren 84 refractive index")
                    elif rf_ice == 1:
                        dir_OP = 'ice_Wrn08/ice_Wrn08'
                        print("Using Warren 08 refractive index")
                    elif rf_ice == 2:
                        dir_OP = 'ice_Pic16/ice_Pic16'
                        print("Using Picard 16 refractive index")

                    FILE_ice = str(dir_mie_ice_files + dir_OP + '_{}.nc'.format(str(grain_rds[i]).rjust(4,'0')))
                    print("\nLayer: {}".format(i))
                    print("Using Mie mode: spheres with radius = {}".format(str(grain_rds[i]).rjust(4,'0')))

            # read in single scattering albedo, MAC and g for ice crystals in each layer,
            # optional with coated liquid water spheres (only available for spherical grains)
            if rwater[i] > grain_rds[i]:

                if grain_shp[i] != 0:
                    raise ValueError("Water coating can only be applied to spherical grains")

                else:
                    # water coating calculations (coated spheres)
                    fn_ice = dir_base + "/Data/rfidx_ice.nc"
                    fn_water = dir_base + "Data/Refractive_Index_Liquid_Water_Segelstein_1981.csv"
                    res = miecoated_driver(rice=grain_rds[i], rwater=rwater[i], fn_ice=fn_ice, rf_ice=rf_ice, fn_water=fn_water, wvl=wvl)
                    SSA_snw[i, :] = res["ssa"]
                    g_snw[i, :] = res["asymmetry"]

                with xr.open_dataset(FILE_ice) as temp:
                    ext_cff_mss = temp['ext_cff_mss'].values
                    MAC_snw[i, :] = ext_cff_mss

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
                    
                    if grain_shp[i] == 0: # snow

                        pass # if layer type is for spheres, no changes required

                    
                    elif grain_shp[i] == 1: # 1 = spheroid, He et al. (2017) parameterization
                        
                        diam_ice = 2.0 * grain_rds[i] # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:

                            fs_sphd = 0.929 # default shape factor for spheroid; He et al. (2017), Table 1
                        
                        else:
                            
                            fs_sphd = shp_fctr[i] # if shp_factor not 0, then use user-defined value
                        
                        
                        if grain_ar[i] == 0:

                            AR_tmp = 0.5 # default aspect ratio for spheroid; He et al. (2017), Table 1
                        
                        else:

                            AR_tmp = grain_ar[i]
                        
                        g_snw_Cg_tmp = g_b0 * (fs_sphd/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * AR_tmp**2 # Eqn. 3.1 in Fu (2007)


                    elif grain_shp[i] == 2: # 3=hexagonal plate, He et al. 2017 parameterization

                        diam_ice = 2.0 * grain_rds[i] # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:
                            
                            fs_hex0 = 0.788 # default shape factor for hexagonal plates; He et al. (2017), Table 1
                        
                        else:

                            fs_hex0 = shp_fctr[i]
                        

                        if grain_ar[i] == 0:

                            AR_tmp = 2.5 # default aspect ratio for hexagonal plate; He et al. (2017), Table 1
                        
                        else:
                        
                            AR_tmp = grain_ar[i]
                                
                        g_snw_Cg_tmp = g_b0 * (fs_hex0/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2   # Eqn. 3.3 in Fu (2007)
                        

                    elif grain_shp[i] == 3: # 4=koch snowflake, He et al. (2017) parameterization

                        diam_ice = 2.0 * grain_rds[i] / 0.544 # effective snow grain diameter
                        
                        if shp_fctr[i] == 0:
                        
                            fs_koch = 0.712 # default shape factor for koch snowflake; He et al. (2017), Table 1
                        
                        else:
                            
                            fs_koch = shp_fctr[i]
                        

                        if grain_ar[i] == 0:

                            AR_tmp = 2.5 # default aspect ratio for koch snowflake; He et al. (2017), Table 1
                        
                        else:

                            AR_tmp = grain_ar[i]
                        
                        g_snw_Cg_tmp = g_b0 * (fs_koch/fs_hex)**g_b1 * diam_ice**g_b2 # Eq.7, He et al. (2017)
                        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 * np.log(AR_tmp) + g_F07_p2 * (np.log(AR_tmp))**2  # Eqn. 3.3 in Fu (2007)
                        
                    
                    if (grain_shp[i] > 0)& (grain_shp[i] < 4):

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
            
            rd = "{}".format(grain_rds[i])
            rd = rd.rjust(4,"0")
            refidx_file = xr.open_dataset(dir_RI_ice+'rfidx_ice.nc')
           
            if rf_ice == 0:
                refidx_re = refidx_file['re_Wrn84'].values
                refidx_im = refidx_file['im_Wrn84'].values 

            elif rf_ice == 1:
                refidx_re = refidx_file['re_Wrn08'].values
                refidx_im = refidx_file['im_Wrn08'].values 

            elif rf_ice == 2:
                refidx_re = refidx_file['re_Pic16'].values
                refidx_im = refidx_file['im_Pic16'].values 

            FILE_ice = str(dir_bubbly_ice + 'bbl_{}.nc').format(rd)
            file = xr.open_dataset(FILE_ice)
            sca_cff_vlm = file['sca_cff_vlm'].values # scattering cross section unit per volume of bubble
            g_snw[i,:] = file['asm_prm'].values
            abs_cff_mss_ice[:] = ((4 * np.pi * refidx_im) / (wvl * 1e-6))/917
            vlm_frac_air = (917 - rho_layers[i]) / 917
            MAC_snw[i,:] = ((sca_cff_vlm * vlm_frac_air) /917) + abs_cff_mss_ice
            SSA_snw[i,:] = ((sca_cff_vlm * vlm_frac_air) /917) / MAC_snw[i,:]

    
    ###################################################
    # Read in impurity optical properties
    ###################################################
    
    # Load optical properties SSA, MAC and g (one row per impurity, one column per wvalengths)
    # Load mass concentrations MSS per layer (one row per layer, one column per umpurity)

    SSAaer = np.zeros([nbr_aer,nbr_wvl])
    MACaer = np.zeros([nbr_aer, nbr_wvl])
    Gaer = np.zeros([nbr_aer,nbr_wvl])
    MSSaer = np.zeros([nbr_lyr, nbr_aer])
    
    for aer in range(nbr_aer):
        impurity_properties = xr.open_dataset(str(dir_mie_lap_files + files[aer]))
        Gaer[aer,:] = impurity_properties['asm_prm'].values
        SSAaer[aer,:] = impurity_properties['ss_alb'].values
        MSSaer[0:nbr_lyr,aer] = mass_concentrations[aer]
        if files[aer] == FILE_brwnC2 or files[aer] == FILE_soot2: #coated particles: use ext_cff_mss_ncl 
            MACaer[aer,:] = impurity_properties['ext_cff_mss_ncl'].values
        else:
            MACaer[aer,:] = impurity_properties['ext_cff_mss'].values


    if GA_units == 1:
        # GA_units ==1 means algae MAC provided in m2/cell
        # this requires conversion into units consistent with other LAPs
        # we assume pigment is the only absorbing material in cell
        # from cell concn we can therefore calculate kg/kg of pigment mix
        
        # MAC in m2/cells/mL: 
        # x 1000 = m2/cells/L
        # / 0.917 = cells/kg  
        # now when mutiplied by n cells, unit is m2/kg
        MACaer[-1,:] = MACaer[-1,:]*1e3/0.917   
        MSSaer[:,:-2] = MSSaer[:,:-2]*1e-9 # skip algae in ng/kg -> kg/kg conversion


    else:  

        MSSaer = MSSaer*1e-9 # all mass concentrations converted to kg/kg unit

    # if the user has provided a C factor, use it to concentrate algae in upper layer
    if isinstance(Cfactor,(int, float)) and (Cfactor > 0):
        MSSaer[0,-1] = MSSaer[0,-1]*Cfactor
        

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

        L_snw[i] = rho_layers[i] * dz[i]
        tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]



    # then for the LAPs in each layer
    for i in range(nbr_lyr):
        for j in range(nbr_aer):

            L_aer[i, j, :] = L_snw[i] * MSSaer[i, j] #kg ice m-2 * cells kg-1 ice = cells m-2
            tau_aer[i, j, :] = L_aer[i, j, :] * MACaer[j, :] # cells m-2 * m2 cells-1

            tau_sum = tau_sum + tau_aer[i, j, :]
            SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
            g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])


    # finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
    for i in range(nbr_lyr):
        
        tau[i,:] = tau_sum[i,:] + tau_snw[i,:]
        SSA[i,:] = (1 / tau[i,:]) * (SSA_sum[i,:] + SSA_snw[i,:] * tau_snw[i,:])
        g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i,:] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))
        
    inputs.tau=tau
    inputs.SSA=SSA
    inputs.g=g
    inputs.L_snw=L_snw
    # just in case any unrealistic values arise (none detected so far)
    SSA[SSA<=0]=0.00000001
    SSA[SSA>=1]=0.99999999
    g[g<=0]=0.00001
    g[g>=1]=0.99999

    # CALL RT SOLVER (TOON  = TOON ET AL, TRIDIAGONAL MATRIX METHOD; 
    # ADD_DOUBLE = ADDING-DOUBLING METHOD)
    
    outputs = c.namedtuple('outputs',['wvl', 'albedo', 'BBA', 'BBAVIS', 'BBANIR', 'abls_slr', 'heat_rt'])

   
    if TOON: 
        
        outputs.wvl, outputs.albedo, outputs.BBA, outputs.BBAVIS, outputs.BBANIR, outputs.abs_slr, outputs.heat_rt = toon_solver(inputs)


    if ADD_DOUBLE:

        outputs.wvl, outputs.albedo, outputs.BBA, outputs.BBAVIS, outputs.BBANIR, outputs.abs_slr, outputs.heat_rt = adding_doubling_solver(inputs)

    return outputs





