from logging import critical


def adding_doubling_solver(inputs):
    

    """
    This script is one of the two optional radiativ transfer solvers available in this package. This
    script deals with the adding-doubling method as translated from MATLAB code from Chloe Whicker 
    (UMich) - October 2020. When it becomes available, any use of this adding-doubling script should cite 
    Chloe's paper.

    This is the appropriate solver for any configuration where solid ice layers and fresnel reflection
    are included.

    """
    
    import numpy as np
    import xarray as xr
    
    #load variables from input table
    tau=inputs.tau
    SSA=inputs.SSA
    g=inputs.g
    nbr_wvl=inputs.nbr_wvl
    wvl=inputs.wvl
    nbr_lyr=inputs.nbr_lyr
    layer_type=inputs.layer_type
    rf_ice=inputs.rf_ice
    R_sfc=inputs.R_sfc
    Fs=inputs.Fs
    Fd=inputs.Fd
    mu_not=inputs.mu_not
    L_snw=inputs.L_snw
    flx_slr=inputs.flx_slr
    dir_base=inputs.dir_base
    DIRECT=inputs.DIRECT
    
    #directory
    dir_RI_ice = str(dir_base + 'Data/') 

    #######################################
    ## DEFINE CONSTANTS AND SET UP ARRAYS
    #######################################
    
    tau0    = tau.T   # read and transpose tau
    g0      = g.T  # read and transpose g
    SSA0  = SSA.T  # read and transpose SSA

    epsilon = 1e-5      # to deal with singularity
    exp_min = 1e-5      # exp(-500)  # minimum number that is not zero - zero will raise error
    trmin   = 1e-5      # minimum transmissivity
    puny    = 1e-10     # not sure how should we define this

    gauspt = [0.9894009, 0.9445750, 0.8656312, 0.7554044, 0.6178762, 0.4580168, 0.2816036, 0.0950125]  # gaussian angles (radians)     
    gauswt = [0.0271525, 0.0622535, 0.0951585, 0.1246290, 0.1495960, 0.1691565, 0.1826034, 0.1894506] # gaussian weights

    vis_max_idx = 50   # index of maximum visible wavelength (0.7 um)
    nir_max_idx = 480 # index of max nir wavelength (5 um)

    # empty arrays
    trndir = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    trntdr = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    trndif = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    rupdir = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    rupdif = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    trndir = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    rdndif = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    fdirup = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    fdifup = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    fdirdn = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    fdifdn = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    dfdir = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    dfdif = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    rdir = np.zeros(shape=[nbr_lyr+1])
    rdif_a = np.zeros(shape=[nbr_lyr+1])
    rdif_b = np.zeros(shape=[nbr_lyr+1])   #layer reflectivity to diffuse radiation from below
    tdir = np.zeros(shape=[nbr_lyr+1])   #layer transmission to direct radiation (solar beam + diffuse)
    tdif_a = np.zeros(shape=[nbr_lyr+1])   #layer transmission to diffuse radiation from above
    tdif_b = np.zeros(shape=[nbr_lyr+1])   #layer transmission to diffuse radiation from below
    trnlay = np.zeros(shape=[nbr_lyr+1])   #solar beam transm for layer (direct beam only)
    F_up = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    F_dwn = np.zeros(shape=[nbr_wvl,nbr_lyr+1])
    F_abs = np.zeros(shape=[nbr_wvl,nbr_lyr])
    F_abs_vis = np.zeros(shape=[nbr_lyr])
    F_abs_nir = np.zeros(shape=[nbr_lyr])
    trndir[:,0] =  1 
    trntdr[:,0] =  1 
    trndif[:,0] =  1 
    rdndif[:,0] =  0 

    # if there are non zeros in layer type, grab the index of the first fresnel layer
    # if there are non-zeros in layer type, load in the precalculated diffuse fresnel reflection
    # (precalculated as large no. of gaussian points required for convergence)
    if np.sum(layer_type) > 0:

        lyrfrsnl = layer_type.index(1)
        print("\nFirst Fresnel bounday is in layer ", lyrfrsnl)

    else:

        lyrfrsnl = 999999999

        # raise error if there are no solid ice layers - in this case use the Toon solver instead!
        print("There are no ice layers in this model configuration\
             - suggest adding a solid ice layer or using faster Toon method")

    # proceed down one layer at a time: if the total transmission to
    # the interface just above a given layer is less than trmin, then no
    # Delta-Eddington computation for that layer is done.
    
    Fresnel_Diffuse_File = xr.open_dataset(dir_RI_ice+'FL_reflection_diffuse.nc')

    for wl in np.arange(0,nbr_wvl,1): # loop through wavelengths
        
        for lyr in np.arange(0,nbr_lyr,1):   # loop through layers

            # open refractive index file and grab real and imaginary parts
            refidx_file = xr.open_dataset(dir_RI_ice+'rfidx_ice.nc')

            if rf_ice == 0:

                refidx_re = refidx_file['re_Wrn84'].values
                refidx_im = refidx_file['im_Wrn84'].values 
                FL_r_dif_a = Fresnel_Diffuse_File['R_dif_fa_ice_Wrn84'].values
                FL_r_dif_b = Fresnel_Diffuse_File['R_dif_fb_ice_Wrn84'].values

            elif rf_ice == 1:
                refidx_re = refidx_file['re_Wrn08'].values
                refidx_im = refidx_file['im_Wrn08'].values 
                FL_r_dif_a = Fresnel_Diffuse_File['R_dif_fa_ice_Wrn08'].values
                FL_r_dif_b = Fresnel_Diffuse_File['R_dif_fb_ice_Wrn08'].values

            elif rf_ice == 2:
                refidx_re = refidx_file['re_Pic16'].values
                refidx_im = refidx_file['im_Pic16'].values
                FL_r_dif_a = Fresnel_Diffuse_File['R_dif_fa_ice_Pic16'].values
                FL_r_dif_b = Fresnel_Diffuse_File['R_dif_fb_ice_Pic16'].values 


            refindx = refidx_re[wl]+refidx_im[wl]  # combine real and imaginary parts into one var
            
            temp1 = refidx_re[wl]**2 - refidx_im[wl]**2 + np.sin(np.arccos(mu_not))**2
            temp2 = refidx_re[wl]**2 - refidx_im[wl]**2 - np.sin(np.arccos(mu_not))**2
            n_real = (np.sqrt(2)/2) * (temp1 + (temp2**2 + 4*refidx_re[wl]**2*refidx_im[wl]**2)**0.5)**0.5

            #  compute next layer Delta-eddington solution only if total transmission
            #  of radiation to the interface just above the layer exceeds trmin.
            
            if trntdr[wl,lyr] > trmin: # condition: only run computation if sufficient
                #flux received from above
                
                mu0 = mu_not  # cosine of beam angle is equal to incident beam
                
                # ice-adjusted real refractive index
                nr = n_real
                # . Eq. 20: Briegleb and Light 2007: adjusts beam angle
                # (i.e. this is Snell's Law for refraction at interface between media)
                # mu0n = -1 represents light travelling vertically upwards and mu0n = +1 
                # represents light travellign vertically downwards
                
                #mu0n = np.sqrt(1-((1-mu0**2)/(refindx*refindx)))  (original, before update for diffuse Fresnel reflection)
                mu0n = np.cos(np.arcsin(np.sin(np.arccos(mu0))/nr))   # this version accounts for diffuse fresnel reflection            
                # condition: if current layer is above fresnel layer or the 
                # top layer is a Fresnel layer
                if lyr < lyrfrsnl or lyrfrsnl==0:
                    
                    mu0n = mu0 

                elif lyr >= lyrfrsnl:
                    # within or below FL
                    mu0n = mu0n 

                # calculation over layers with penetrating radiation
                # includes optical thickness, single scattering albedo, 
                # asymmetry parameter and total flux
                tautot = tau0[wl,lyr] 
                wtot   = SSA0[wl,lyr] 
                gtot   = g0[wl,lyr] 
                ftot   = g0[wl,lyr] * g0[wl,lyr] 
                
                # coefficient for delta eddington solution for all layers 
                # Eq. 50: Briegleb and Light 2007
                ts   = (1-(wtot * ftot)) * tautot # layer delta-scaled extinction optical depth
                ws   = ((1-ftot) * wtot) / (1-(wtot * ftot)) # layer delta-scaled single scattering albedo
                gs   = (gtot-ftot)/(1-ftot) # layer delta-scaled asymmetry parameter
                lm   = np.sqrt(3 * (1-ws) * (1-ws * gs)) # lambda
                ue   = 1.5 * (1-ws * gs) / lm # u equation, term in diffuse reflectivity and transmissivity
                
                extins = max(exp_min, np.exp(-lm * ts)) # extinction, MAX function lyr keeps from getting an error if the exp(-lm*ts) is < 1e-5
                ne = (ue+1)**2 / extins - (ue-1)**2 * extins # N equation, term in diffuse reflectivity and transmissivity
                
                 # ! first calculation of rdif, tdif using Delta-Eddington formulas
                 # Eq.: Briegleb 1992  alpha and gamma for direct radiation

                rdif_a[lyr] = (ue**2-1) * (1/extins - extins)/ne # R BAR = layer reflectivity to DIFFUSE radiation
                tdif_a[lyr] = 4*ue/ne # T BAR layer transmissivity to DIFFUSE radiation
                
                 # evaluate rdir, tdir for direct beam
                trnlay[lyr] = max([exp_min, np.exp(-ts/mu0n)]) # transmission from TOA to interface
                
                 #  Eq. 50: Briegleb and Light 2007  alpha and gamma for direct radiation
                alp = (0.75 * ws * mu0n) * ((1 + gs * (1-ws)) / (1 - lm**2 * mu0n**2 + epsilon))   #alp = alpha(ws,mu0n,gs,lm)
                gam = (0.5 * ws) * ((1 + 3 * gs * mu0n**2 * (1-ws)) / (1-lm**2 * mu0n**2 + epsilon))     #gam = gamma(ws,mu0n,gs,lm)
                
                # apg = alpha plus gamma
                # amg = alpha minus gamma
                apg = alp + gam 
                amg = alp - gam 

                rdir[lyr] = apg*rdif_a[lyr] +  amg*(tdif_a[lyr]*trnlay[lyr] - 1)     #layer reflectivity to DIRECT radiation
                tdir[lyr] = apg*tdif_a[lyr] + (amg* rdif_a[lyr]-apg+1)*trnlay[lyr]   #layer transmissivity to DIRECT radiation
                
                 # recalculate rdif,tdif using direct angular integration over rdir,tdir,
                 # since Delta-Eddington rdif formula is not well-behaved (it is usually
                 # biased low and can even be negative)  use ngmax angles and gaussian
                 # integration for most accuracy:
                
                R1 = rdif_a[lyr]   # use R1 as temporary var
                T1 = tdif_a[lyr]   # use T1 as temporary var
                swt = 0 
                smr = 0 
                smt = 0 
                
                 # loop through the gaussian angles for the AD integral
                for ng in np.arange(0,len(gauspt),1):     #gaussian angles (radians)
                    
                    mu  = gauspt[ng]         # solar zenith angles
                    gwt = gauswt[ng]         # gaussian weight
                    swt = swt + mu*gwt       # sum of weights
                    trn = max([exp_min, np.exp(-ts/mu)])   # transmission
                    
                    alp = (0.75*ws*mu) * (1 + gs * (1-ws)) / (1 - lm**2 * mu**2 + epsilon)   #alp = alpha(ws,mu0n,gs,lm)
                    gam = (0.5 * ws) * (1 + 3 * gs * mu**2 * (1-ws)) / (1-lm**2 * mu**2 + epsilon)  #gam = gamma(ws,mu0n,gs,lm)
                    
                    apg = alp + gam 
                    amg = alp - gam 
                    rdr = apg*R1 + amg*T1*trn - amg 
                    tdr = apg*T1 + amg*R1*trn - apg*trn + trn 
                    smr = smr + mu*rdr*gwt   #accumulator for rdif gaussian integration
                    smt = smt + mu*tdr*gwt   #accumulator for tdif gaussian integration

                
                rdif_a[lyr] = smr/swt 
                tdif_a[lyr] = smt/swt 
                
                #! homogeneous layer
                rdif_b[lyr] = rdif_a[lyr] 
                tdif_b[lyr] = tdif_a[lyr] 
                
                ###################################################
                # Fresnel layer
                ##############################################################
                
                if lyr == lyrfrsnl:
                    
                    refindx = np.complex(refidx_re[wl],refidx_im[wl])
                    critical_angle = np.arcsin(refindx)


                    if np.arccos(mu_not) < critical_angle:
                        # in this case, no total internal reflection

                        #! compute fresnel reflection and transmission amplitudes
                        #! for two polarizations: 1=perpendicular and 2=parallel to
                        #! the plane containing incident, reflected and refracted rays.
                        
                        #! Eq. 22  Briegleb & Light 2007
                        # inputs to equation 21 (i.e. Fresnel formulae for R and T)
                        R1 = (mu0-nr*mu0n) / (mu0 + nr*mu0n)    #reflection amplitude factor for perpendicular polarization
                        R2 = (nr*mu0 - mu0n) / (nr*mu0 + mu0n)    #reflection amplitude factor for parallel polarization
                        T1 = 2*mu0 / (mu0 + nr*mu0n)                   #transmission amplitude factor for perpendicular polarization
                        T2 = 2 * mu0 / (nr * mu0 + mu0n)                   #transmission amplitude factor for parallel polarization
                        
                        #! unpolarized light for direct beam
                        #! Eq. 21  Brigleb and light 2007
                        Rf_dir_a = 0.5 * (R1**2 + R2**2) 
                        Tf_dir_a = 0.5 * (T1**2 + T2**2) * nr * mu0n / mu0 
                    
                    else: # in this case, total internal reflection occurs
                        Tf_dir_a = 0
                        Rf_dir_a = 1
                    
                    # precalculated diffuse reflectivities and transmissivities
                    # for incident radiation above and below fresnel layer, using
                    # the direct albedos and accounting for complete internal
                    # reflection from below. Precalculated because high order
                    # number of gaussian points (~256) is required for convergence:
                    
                    # Eq. 25  Brigleb and light 2007
                    # diffuse reflection of flux arriving from above
                    
                    Rf_dif_a = FL_r_dif_a[wl]             # reflection from diffuse unpolarized radiation
                    Tf_dif_a = 1 - Rf_dif_a     #t ransmission from diffuse unpolarized radiation
                    
                    # diffuse reflection of flux arriving from below
                    Rf_dif_b = FL_r_dif_b[wl] 
                    Tf_dif_b = 1 - Rf_dif_b 
                    
                    ######################################################################
                    # the lyr = lyrfrsnl layer properties are updated to combine
                    # the fresnel (refractive) layer, always taken to be above
                    # the present layer lyr (i.e. be the top interface):
                    
                    rintfc   = 1 / (1-Rf_dif_b*rdif_a[lyr])  # denom interface scattering
                    
                    # layer transmissivity to DIRECT radiation
                    # Eq. B7  Briegleb & Light 2007
                    tdir[lyr] = Tf_dir_a * tdir[lyr] + Tf_dir_a*rdir[lyr] * Rf_dif_b*rintfc*tdif_a[lyr]   

                    # layer reflectivity to DIRECT radiation
                    # Eq. B7  Briegleb & Light 2007
                    rdir[lyr] = Rf_dir_a + Tf_dir_a*rdir[lyr] * rintfc * Tf_dif_b 
                    
                    # R BAR = layer reflectivity to DIFFUSE radiation (above)
                    # Eq. B9  Briegleb & Light 2007
                    rdif_a[lyr] = Rf_dif_a + Tf_dif_a*rdif_a[lyr] * rintfc * Tf_dif_b 
                    
                    # R BAR = layer reflectivity to DIFFUSE radiation (below)
                    # Eq. B10  Briegleb & Light 2007
                    rdif_b[lyr] = rdif_b[lyr] + tdif_b[lyr] * Rf_dif_b * rintfc * tdif_a[lyr] 
                    
                    # T BAR layer transmissivity to DIFFUSE radiation (above),
                    # Eq. B9  Briegleb & Light 2007
                    tdif_a[lyr] = tdif_a[lyr] * rintfc * Tf_dif_a   

                    # Eq. B10  Briegleb & Light 2007
                    tdif_b[lyr] = tdif_b[lyr] * rintfc * Tf_dif_b   
                    
                    #! update trnlay to include fresnel transmission
                    trnlay[lyr] = Tf_dir_a*trnlay[lyr] 
                    
                    # end lyr = lyrfrsnl condition
                    # end trntdr[lyr, wl] > trmin condition
                
                #  ! Calculate the solar beam transmission, total transmission, and
                #  ! reflectivity for diffuse radiation from below at interface lyr,
                #  ! the top of the current layer lyr:
                #  !
                #  !              layers       interface
                #  !
                #  !       ---------------------  lyr-1
                #  !                lyr-1
                #  !       ---------------------  lyr
                #  !                 lyr
                #  !       ---------------------
                #  ! note that we ignore refraction between sea ice and underlying ocean:
                #  !
                #  !              layers       interface
                #  !
                #  !       ---------------------  lyr-1
                #  !                lyr-1
                #  !       ---------------------  lyr
                #  !       \\\\\\\ ocean \\\\\\\
                
            
            # Eq. 51  Briegleb and Light 2007

            trndir[wl,lyr+1] = trndir[wl,lyr]*trnlay[lyr]  # solar beam transmission from top
            # trnlay = exp(-ts/mu_not) = direct solar beam transmission

            # interface multiple scattering for lyr-1
            refkm1 = 1/(1 - rdndif[wl,lyr]*rdif_a[lyr])  
            
            # direct tran times layer direct ref
            tdrrdir = trndir[wl,lyr]*rdir[lyr]            
            
            # total down diffuse = tot tran - direct tran
            tdndif = trntdr[wl,lyr] - trndir[wl,lyr]     
            
            # total transmission to direct beam for layers above
            trntdr[wl,lyr+1] = trndir[wl,lyr]*tdir[lyr] + (tdndif + tdrrdir*rdndif[wl,lyr])*refkm1*tdif_a[lyr] 
            
            # Eq. B4  Briegleb and Light 2007
            rdndif[wl,lyr+1] = rdif_b[lyr] + (tdif_b[lyr]*rdndif[wl,lyr]*refkm1*tdif_a[lyr])    #reflectivity to diffuse radiation for layers above
            trndif[wl,lyr+1] = trndif[wl,lyr]*refkm1*tdif_a[lyr]   #diffuse transmission to diffuse beam for layers above
        
        # end main level loop  number of layers
        
        # ! compute reflectivity to direct and diffuse radiation for layers
        # ! below by adding succesive layers starting from the underlying
        # ! ocean and working upwards:
        # !
        # !              layers       interface
        # !
        # !       ---------------------  lyr
        # !                 lyr
        # !       ---------------------  lyr+1
        # !                lyr+1
        # !       ---------------------
        
        # set the underlying ground albedo
        rupdir[wl,nbr_lyr] = R_sfc[wl]    # reflectivity to direct radiation for layers below
        rupdif[wl,nbr_lyr] = R_sfc[wl]    # reflectivity to diffuse radiation for layers below


        for lyr in np.arange(nbr_lyr-1,-1,-1):  # starts at the bottom and works its way up to the top layer
            
            #Eq. B5  Briegleb and Light 2007
            #! interface scattering
            refkp1 = 1/( 1 - rdif_b[lyr]*rupdif[wl,lyr+1]) 
            
            # dir from top layer plus exp tran ref from lower layer, interface
            # scattered and tran thru top layer from below, plus diff tran ref
            # from lower layer with interface scattering tran thru top from below
            rupdir[wl,lyr] = rdir[lyr] + (trnlay[lyr] * rupdir[wl,lyr+1] + (tdir[lyr]-trnlay[lyr])* rupdif[wl,lyr+1])*refkp1*tdif_b[lyr] 
            
            # dif from top layer from above, plus dif tran upwards reflected and
            # interface scattered which tran top from below
            rupdif[wl,lyr] = rdif_a[lyr] + tdif_a[lyr]*rupdif[wl,lyr+1]*refkp1*tdif_b[lyr] 

    
    # fluxes at interface

    for wl in np.arange(0,nbr_wvl,1):
        
        for lyr in np.arange(0,nbr_lyr+1,1):
            
            # Eq. 52  Briegleb and Light 2007
            # interface scattering
            refk = 1/(1 - rdndif[wl,lyr]*rupdif[wl,lyr]) 
            
            # dir tran ref from below times interface scattering, plus diff
            # tran and ref from below times interface scattering
            fdirup[wl,lyr] = (trndir[wl,lyr]*rupdir[wl,lyr] + (trntdr[wl,lyr]-trndir[wl,lyr]) * rupdif[wl,lyr])*refk 
            
            # dir tran plus total diff trans times interface scattering plus
            # dir tran with up dir ref and down dif ref times interface scattering
            fdirdn[wl,lyr] = trndir[wl,lyr] + (trntdr[wl,lyr]- trndir[wl,lyr] + trndir[wl,lyr] * rupdir[wl,lyr] * rdndif[wl,lyr])*refk 
            
            # diffuse tran ref from below times interface scattering
            fdifup[wl,lyr] = trndif[wl,lyr]*rupdif[wl,lyr]*refk 
            
            # diffuse tran times interface scattering
            fdifdn[wl,lyr] = trndif[wl,lyr]*refk 
            
            # dfdir = fdirdn - fdirup
            dfdir[wl,lyr] = trndir[wl,lyr] + (trntdr[wl,lyr]-trndir[wl,lyr]) * (1 - rupdif[wl,lyr]) * refk - trndir[wl,lyr]*rupdir[wl,lyr]\
                * (1 - rdndif[wl,lyr]) * refk 
            
            if dfdir[wl,lyr] < puny:
            
                dfdir[wl,lyr] = 0  #!echmod necessary?
                # dfdif = fdifdn - fdifup
            
            dfdif[wl,lyr] = trndif[wl,lyr] * (1 - rupdif[wl,lyr]) * refk 
            
            if dfdif[wl,lyr] < puny:

                dfdif[wl,lyr] = 0  #!echmod necessary?


    # ----- End Radiative Solver Adding Doubling Method -----
    # ----- Calculate fluxes ----


    for n in np.arange(0,nbr_lyr+1,1):
        
        F_up[:,n]  = (fdirup[:,n]*(Fs*mu_not*np.pi) + fdifup[:,n]*Fd) 
        F_dwn[:,n] = (fdirdn[:,n]*(Fs*mu_not*np.pi) + fdifdn[:,n]*Fd) 

    F_net = F_up - F_dwn 

    # Absorbed flux in each layer
    #F_abs[:,0:nbr_lyr] = F_net[:,1:]-F_net[:,0:-1]
    F_abs[:,:] = F_net[:,1:]-F_net[:,0:-1]

    # albedo
    acal  = F_up[:,0]/F_dwn[:,0] 

    # Upward flux at upper model boundary
    F_top_pls = F_up[:,0] 
    
    # Net flux at lower model boundary = bulk transmission through entire
    # media = absorbed radiation by underlying surface:
    F_btm_net = -F_net[:,nbr_lyr] 

    # Spectrally-integrated absorption in each layer:
    F_abs_slr = np.sum(F_abs,axis=0) 
    
    for i in np.arange(0,nbr_lyr,1):

        F_abs_vis[i] = sum(F_abs[0:vis_max_idx,i]) 
        F_abs_nir[i] = sum(F_abs[vis_max_idx:nir_max_idx,i]) 

    # Spectrally-integrated absorption by underlying surface:
    F_abs_btm = np.sum(F_btm_net,axis=0) 
    F_abs_vis_btm = np.sum(F_btm_net[0:vis_max_idx],axis=0) 
    F_abs_nir_btm = np.sum(F_btm_net[vis_max_idx:nir_max_idx+1],axis=0) 

    # Radiative heating rate:
    heat_rt = F_abs_slr/(L_snw*2117)    #[K/s] 2117 = specific heat ice (J kg-1 K-1)
    heat_rt = heat_rt*3600               #[K/hr]

    # Energy conservation check:
    # Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
    energy_sum = (mu_not*np.pi*Fs)+Fd - (np.sum(F_abs,axis=1) + F_btm_net + F_top_pls) 

    energy_conservation_error = sum(abs(energy_sum)) 

    if energy_conservation_error > 1e-10:

        print('energy conservation error: {}'.format(energy_conservation_error))

    # Hemispheric wavelength-dependent albedo:
    if DIRECT ==1: 
    
        albedo = F_top_pls/((mu_not*np.pi*Fs)+Fd) 
    
    else:
        albedo = rupdif[:,0] 


    #double check if the albedo calculated are the same
    adif = np.sum(acal - albedo) 
    
    if adif > 1e-10:
        
        print('error in albedo calculation')

    albedo = acal 


    # Spectrally-integrated solar, visible, and NIR albedos:
        
    alb_bb = np.sum(flx_slr*albedo)/np.sum(flx_slr) 

    alb_vis = np.sum(flx_slr[0:vis_max_idx] * albedo[0:vis_max_idx]) / np.sum(flx_slr[0:vis_max_idx]) 

    alb_nir = np.sum(flx_slr[vis_max_idx:nir_max_idx] * albedo[vis_max_idx:nir_max_idx]) / np.sum(flx_slr[vis_max_idx:nir_max_idx]) 

    # Spectrally-integrated VIS and NIR total snowpack absorption:
    abs_vis = np.sum(flx_slr[0:vis_max_idx] * (1-albedo[0:vis_max_idx])) 
    abs_nir = np.sum(flx_slr[vis_max_idx:nir_max_idx]*(1-albedo[vis_max_idx:nir_max_idx])) 

    #########################  OUTPUT  #############################

    flx_dwn_spc = mu_not*np.pi*Fs+Fd  # spectral downwelling flux at model top [W/m2/band]
    alb_slr = alb_bb              # solar broadband albedo
    abs_snw_slr = np.sum(F_abs_slr)      # total solar absorption by entire snow column (not including underlying substrate) [W/m2]
    # abs_snw_vis = np.sum(F_abs_vis)      # visible solar absorption by entire snow column (not including underlying substrate) [W/m2]
    # abs_snw_nir = np.sum(F_abs_nir)      # near-IR solar absorption by entire snow column (not including underlying substrate) [W/m2]
    # abs_spc = np.sum(F_abs,axis=1)      # spectral absorption by entire snow column [W/m2/band]

    # abs_snw_top_slr = F_abs_slr[0]        # top snow layer solar absorption [W/m2]
    # abs_snw_top_vis = F_abs_vis[0]        # top snow layer VIS absorption [W/m2]
    # abs_snw_top_nir = F_abs_nir[0]        # top snow layer NIR absorption [W/m2]

    # abs_ground_slr  = F_abs_btm           # total solar absorption by underlying substrate [W/m2]
    # abs_ground_vis  = F_abs_vis_btm       # visible absorption by underlying substrate [W/m2]
    # abs_ground_nir  = F_abs_nir_btm       # near-IR absorption by underlying substrate [W/m2]


    return wvl, albedo, alb_bb, alb_vis, alb_nir, F_abs_slr, heat_rt
