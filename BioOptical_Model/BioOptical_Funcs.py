#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Contains functions relating to bio-optical model components of BioSNICAR_GO_py
Run from BioOptical_driver.py

"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from scipy.signal import savgol_filter
import scipy as sci
from scipy.integrate import quad
import csv
from math import exp
from miepython import mie



def bio_optical(wdir = '/home/joe/Code/BioSNICAR_GO_PY/', load_MAC = True, apply_packaging_correction=True, calc_MAC = False, calc_k = True, pig_mass = True,
                pig_frac = False, Pottier = False, Cook = True, cell_dm_weight = 0.82, chla = 0.01, chlb = 0.00066,
                ppro = 0.01, psyn = 0, purp = 0.068, Xw = 0.8, density= 1400, nm = 1.4, smooth = True, smoothStart = 44,
                smoothStop = 100, window_size = 25, poly_order = 3, savefiles = False, savepath = "path", savefilename = "name", 
                plot_optical_props = True, plot_pigment_MACs = True, saveplots = True):


    ##########
    # SET UP
    ##########
    
    # set up empty containers
    data = pd.DataFrame() # set up dataframe for storing optical properties
    WL = []
    Ea1 = []
    Ea2 = []
    Ea3 = []
    Ea4 = []
    Ea5 = []
    Ea1n = []
    Ea2n = []
    Ea3n = []
    Ea4n = []
    Ea5n = []
    WatRI = []
    WatRIn = []
    k_list = []
    real_list = []
    MACcorr = []
    MACcorrn = []

    #################################
    # CASE 1: MAC is loaded from file
    #################################

    if load_MAC: # choose this option to load an empirically derived MAC from file

        if apply_packaging_correction: # optionally apply correction for packaging effect
            # note that this option was included only for experiments related to the 
            # optical effects of the phenol pigment in isolation for Williamson et al (2020).
            # Usually, an AMC loaded from file will be for an entire cell, so no packaging
            # effect required.
            
            # read in preprocessed file containing MAC corrected for packaging effects
            MAC = pd.read_csv(str(wdir+'Data/phenol_mac_packaging_corrected.csv'),header=None,names=['MAC'])
            MAC = MAC[0:4695] # subsample to appropriate resolution for snicar
            MAC = MAC[0:-1:10] 
            data['MAC'] = MAC['MAC'].dropna() # drop NaNs and save to dataframe

        else:
            # this is the more likely option for loaded MACs - the MAC likely represents the 
            # total cell MAC measured empirically
            MAC = pd.read_csv(str(wdir+'Data/Empirical_MAC.csv'),header=None,names=['MAC'])
            MAC = MAC[0:4695] # subsample to appropriate resolution for snicar
            MAC = MAC[0:-1:10]
            data['MAC'] = MAC['MAC'].dropna() # drop NaNs and save to dataframe

    
    ##############################################################
    # CASE 2: MAC and k are calculated from pigment concentrations
    ##############################################################
    # LOAD PIGMENT DATA
    ###################

    if calc_MAC or calc_k: # choose this option to calculate MAC and/or k theoretically
        
        # set wavelengths
        WL = np.arange(300,5000,10) # in nanometers
        WLmeters = [i*1e-9 for i in (WL)] # in meters

        # open in-vivo absorption coefficients (m2/mg) for each algal pigment

        # 1) Chlorophyll a
        with open(str(wdir+'Data/chlorophyll-a.csv'))as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea1.append(cellf)
        
        # 2) Chlorophyll b
        with open(str(wdir+'Data/chlorophyll-b.csv')) as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea2.append(cellf)
        
        # 3) Photoprotective carotenoids
        with open(str(wdir+'Data/Photoprotective_carotenoids.csv'))as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea3.append(cellf)

        # 4) Photosynthetic carotenoids
        with open(str(wdir+'Data/Photosynthetic_carotenoids.csv'))as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea4.append(cellf)

        # 5) Phenol (there are multiple files for the phenol available in
        # the repository representing different measurement conditions.
        # Recommend using "inVivoPhenolMAC.csv as default")
        with open(str(wdir+'Data/inVivoPhenolMAC.csv'))as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea5.append(cellf)

        # Water
        with open(str(wdir+'Data/water_RI.csv'))as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    WatRI.append(cellf)
      
        # In case user wishes to apply packaging correction, load the 
        # vector of correction factors to be applied to the MACs at each
        # wavelength (see Williamson et al 2020)
        with open(str(wdir+'Data/phenol_mac_correction.csv'))as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    MACcorr.append(cellf)


        ########################
        # Pigment PreProcessing
        ########################

        # ensure all pigment MACs are equal in length and cover identical
        # wavelength range, and that their spectral resolution is appropriate
        # for integrating into the RT model.

        for i in np.arange(len(MACcorr),750,1): # ensure list length = 750
            MACcorr.append(1)

        # extend pigment data down to 300 nm by padding with the value at 350nm and
        # small but non-zero value from 750 - 5000 nm (avoid divide by zero errors)    
        Ea1 = [Ea1[0] for _ in range(50)] + Ea1
        Ea2 = [Ea2[0] for _ in range(50)] + Ea2
        Ea3 = [Ea3[0] for _ in range(50)] + Ea3
        Ea4 = [Ea4[0] for _ in range(50)] + Ea4
        
        for i in np.arange(751,5000,1):
            Ea1.append(0)
            Ea2.append(0)
            Ea3.append(0)
            Ea4.append(0)
            Ea5.append(0)
            MACcorr.append(1)
    
        # downsample pigment data to match SNICAR resolution) 
        for i in np.arange(1,len(Ea1),10):
            Ea1n.append(Ea1[i])
            Ea2n.append(Ea2[i])
            Ea3n.append(Ea3[i])
            Ea4n.append(Ea4[i])
            Ea5n.append(Ea5[i])
            WatRIn.append(WatRI[i])
            MACcorrn.append(MACcorr[i])

        if apply_packaging_correction:
            Ea5n = [Ea5n[i] * MACcorrn[i] for i in np.arange(0,len(MACcorrn)-1,1)]

        # copy water RI over zeros at non-absorbing wavelengths so that 
        # cells look like water at non-absorbing wavelengths
        Ea1n[44:-1] = WatRIn[44:-1]
        Ea2n[44:-1] = WatRIn[44:-1]
        Ea3n[44:-1] = WatRIn[44:-1]
        Ea4n[44:-1] = WatRIn[44:-1]
        Ea5n[44:-1] = WatRIn[44:-1]


    ################
    # Calculate MAC
    ################

    if calc_MAC: # if MAC is calculated on-the-fly using theoretical approach
        
        if pig_frac:
        # if data was provided as mass fraction convert percentage dw pigment to actual mass of pigment per cell (mg)
            chla_w = chla * cell_dm_weight*1e-6
            chlb_w = chlb * cell_dm_weight*1e-6
            ppro_w = ppro * cell_dm_weight*1e-6
            psyn_w = psyn * cell_dm_weight*1e-6
            purp_w = purp * cell_dm_weight*1e-6

        elif pig_mass:
        # If data was provided in units of mg pigment/cell, read directly from user values.
            chla_w = chla
            chlb_w = chlb
            ppro_w = ppro
            psyn_w = psyn
            purp_w = purp

        # Multiply mass of each pigment (mg) by absorption coefficient (m2/mg)
        # at each wavelenth giving absorption in units of m2 per cell
        EW1m = [a * chla_w for a in Ea1n]
        EW2m = [a * chlb_w for a in Ea2n]
        EW3m = [a * ppro_w for a in Ea3n]
        EW4m = [a * psyn_w for a in Ea4n]
        EW5m = [a * purp_w for a in Ea5n]

        EWWm = [sum(x) for x in zip(EW1m,EW2m,EW3m,EW4m,EW5m)] # Sum all pigments (m2)
        MAC = [(i/cell_dm_weight)*1e12 for i in EWWm] # normalise from units of cells to units of mass (x 1e12 to get m2/kg)
        data['MAC'] = MAC # save MAC to dataframe

    ################################
    # Calculate imaginary part of RI
    ################################

    if calc_k: # to calculate imaginary refractive index for the cell
        
        if pig_frac: # if pigment values are provided as a mass fraction
        # follow Pottier 2005 / Dauchet (2015) / Cook et al (2017) route to imaginary RI
        
        # multiply pigment MAC by % dw and convert from m2/mg to m2/kg
           EW1 = [a * 1e6 * chla for a in Ea1n]
           EW2 = [a * 1e6 * chlb for a in Ea2n]
           EW3 = [a * 1e6 * ppro for a in Ea3n]
           EW4 = [a * 1e6 * psyn for a in Ea4n]
           EW5 = [a * 1e6 * purp for a in Ea5n]
           EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)] # Sum all pigments       

        else: # if pigment values are provided as absolute mass (mg)
        # if data provided in mg pigment per cell, divide by cell weight in mg
        # to get mass fraction. multiply by pigment MAC (x 1e6 = m2/kg)

            chla_frac = chla/(cell_dm_weight*1e-6)
            chlb_frac = chlb/(cell_dm_weight*1e-6)
            ppro_frac = ppro/(cell_dm_weight*1e-6)
            psyn_frac = psyn/(cell_dm_weight*1e-6)
            purp_frac = purp/(cell_dm_weight*1e-6)

            EW1 = [a  * 1e6 * chla_frac for a in Ea1n]
            EW2 = [a  * 1e6 * chlb_frac for a in Ea2n]
            EW3 = [a  * 1e6 * ppro_frac for a in Ea3n]
            EW4 = [a  * 1e6 * psyn_frac for a in Ea4n]
            EW5 = [a  * 1e6 * purp_frac for a in Ea5n]
            EWW = [sum(x) for x in zip(EW1,EW2,EW3,EW4,EW5)] # Sum all pigments

        # Calculate imaginary refractive index (k)
        # according to user choice of Pottier (2005) method or Cook (2017) method
        
        for i in np.arange(0,len(WL),1):
        
            if Pottier: # Pottier et al (2005) equation
                k = (((1 - Xw) / Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) #original Pottier equation

            elif Cook: # Cook et al (2019) update
                k = (Xw * WatRIn[i]) + ((1 - Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) # Cook (2017) updated version

            k_list.append(k)
            real_list.append(nm)

        data['Imag'] = k_list # append imaginary RI to dataframe
        data['Real'] = real_list # append real RI to dataframe
    

    ######################################
    # Post-processing, saving and plotting
    ######################################

    # apply optional smoothing filter with user defined window width and polynomial order
    if smooth:
        yhat = savgol_filter(MAC[smoothStart:smoothStop], window_size, poly_order) # window size 51, polynomial order 3
        MAC[smoothStart:smoothStop] = yhat

    # optionally save files to savepath
    if savefiles: # optional save dataframe to csv files
        data.to_csv(str(savepath+'/Data/'+'{}_Dataset.csv'.format(savefilename)))
        data['Imag'].to_csv(str(savepath+'/Data/'+'{}KK.csv'.format(savefilename)),header=None,index=False)
        data['MAC'].to_csv(str(savepath+'/Data/'+'{}MAC.csv'.format(savefilename)),header=None,index=False)
        data['Real'].to_csv(str(savepath+'/Data/'+'{}Real.csv'.format(savefilename)),header=None,index=False)

    # optionally plot figures to interative window
    if plot_optical_props:
        plt.figure(figsize=(10,15))
        plt.subplot(2,1,1)
        plt.plot(WL[0:220],MAC[0:220]),plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('MAC (m^2/kg)',fontsize=16)
        plt.tight_layout()

        plt.subplot(2,1,2)
        plt.plot(WL[0:220],k_list[0:220]), plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('K (dimensionless)',fontsize=16)
        plt.tight_layout()

        if saveplots:
            plt.show()
            plt.savefig(str(savepath+"/MACandK.png"))

        else:
            plt.show()
        
    if plot_pigment_MACs:

        plt.figure(figsize=(12, 12))
        plt.plot(WL, Ea1n, label='Chlorophyll a')
        plt.plot(WL, Ea2n, label='Chlorophyll b')
        plt.plot(WL, Ea3n, label='Secpndary carotenoids')
        plt.plot(WL, Ea4n, label='Primary carotenoids')
        plt.plot(WL, Ea5n, label='Purpurogallin')
        plt.xlabel('Wavelengths nm')
        plt.ylabel('In vivo mass absorption coefficient (m2/mg)')
        plt.xlim(300, 750)

        plt.legend(loc='best')

        if saveplots:
            plt.show()
            plt.savefig(str(savepath+"/PigmentMACs.png"))
        
        else:
            plt.show()


    return k_list, real_list, MAC, data


def preprocess_RI(RealsPath,ImagsPath,MACPath):
    """
    Function for preprocessing inpout data to geometrical optics calculations

    parameters:
    RealsPath: path to real part of RI
    ImagsPath: path to imaginary part of RI
    MACPath: path to mass absorption coefficient

    returns:
    reals: numpy array of preprocessed real part of RI
    imags: numpy array of preprocessed imaginary part of RI
    MAC: numpy array of preprocessed MAC

    """
    # load input files for real RI, imaginary RI, mass absorption coefficient. Set wavelength array.
    wavelengths = np.arange(0.305, 5, 0.01)  # 300 - 5000 nm in 10 nm steps
    reals = pd.read_csv(RealsPath, header=None)
    reals = np.array(reals)
    imags = pd.read_csv(ImagsPath, header=None)
    imags = np.array(imags)
    MAC = pd.read_csv(MACPath, names=['vals'], header=None, index_col=None)
    MAC = np.array(MAC['vals'])

    return reals, imags, MAC, wavelengths


def calc_optical_params_GO(savepath, r, depth, reals, imags, wavelengths, plots=False, savefigs = True, report_dims=False):

    """
    Function calculates single scattering optical properties using geometric optics approximation of van Diederhovn (2014)
    assuming a hexagonal cylindrical shape.

    params:
    savepath: path for saving plots
    r: radius of cell cross section
    depth: length of cell along long axis
    reals: numpy array of preprocessed real part of RI
    imags: numpy array of preprocessed imaginary part of RI
    wavelengths: numpy array of wavelengths at appropriate range and resolution
    plots: Boolean toggling plotting on/off
    savefigs: Boolean toggling figure saving on/off
    report_dims: Boolean toggling cell dimensions printing to console

    returns:

    Assy_list: list of assymetry parameter per unit wavelength
    SSA_list: list of single scattering albedo per unit wavelength
    absXS_list: list of absorption cross section per unit wavelength
    MAC: list of mass absorption coefficient per unit wavelength
    depth: scalar for cell length
    r: scalar for cell radius
    Chi_abs_list: list of absorption coefficient per unit wavelength
    Reff: effective radius (i.e. radius of sphere with equal SA:vol as cylindrical cell)
    X_list: list of size parameter X
    
    """
    
    # set up lists
    SSA_list = []
    Assy_list = []
    absXS_list = []
    Chi_abs_list = []
    X_list = []

    # calculate cylinder dimensions
    diameter = 2 * r
    V = depth * (np.pi * r ** 2)  # volume of cylinder
    Reff = (V / ((
                             4 / 3) * np.pi)) ** 1 / 3  # effective radius (i.e. radius of sphere with equal volume to real cylinder)
    Area_total = 2 * (np.pi * r ** 2) + (2 * np.pi * r) * (
        depth)  # total surface area - 2 x ends plus circumference * length
    Area = Area_total / 4  # projected area
    ar = diameter / depth  # aspect ratio
    delta = 0.3

    # start geometric optics calculation per wavelength
    for i in np.arange(0, len(wavelengths), 1):
        mr = reals[i]  # load real RI per wavelength
        mi = imags[i]  # load imaginary RI per wavelength
        wl = wavelengths[i]  # loadd wavelength

        # SSA parameterization
        a = [0.457593, 20.9738]  # for ar=1

        # SSA correction for AR != 1 (Table 2)
        nc1 = 3
        nc2 = 4
        c_ij = np.zeros(nc1 * nc2 * 2).reshape((nc1, nc2, 2))
        #   ---------- Plates ----------
        c_ij[:, 0, 0] = [0.000527060, 0.309748, -2.58028]
        c_ij[:, 1, 0] = [0.00867596, -0.650188, -1.34949]
        c_ij[:, 2, 0] = [0.0382627, -0.198214, -0.674495]
        c_ij[:, 3, 0] = [0.0108558, -0.0356019, -0.141318]
        #   --------- Columns ----------
        c_ij[:, 0, 1] = [0.000125752, 0.387729, -2.38400]
        c_ij[:, 1, 1] = [0.00797282, 0.456133, 1.29446]
        c_ij[:, 2, 1] = [0.00122800, -0.137621, -1.05868]
        c_ij[:, 3, 1] = [0.000212673, 0.0364655, 0.339646]

        # diffraction g parameterization
        b_gdiffr = [-0.822315, -1.20125, 0.996653]

        # raytracing g parameterization ar=1
        p_a_eq_1 = [0.780550, 0.00510997, -0.0878268, 0.111549, -0.282453]

        # ---- g correction for AR != 1 (Also applied to AR=1 as plate) (Table 3)
        nq1 = 3
        nq2 = 7
        q_ij = np.zeros(nq1 * nq2 * 2).reshape((nq1, nq2, 2))
        #   ---------- Plates ----------
        q_ij[:, 0, 0] = [-0.00133106, -0.000782076, 0.00205422]
        q_ij[:, 1, 0] = [0.0408343, -0.00162734, 0.0240927]
        q_ij[:, 2, 0] = [0.525289, 0.418336, -0.818352]
        q_ij[:, 3, 0] = [0.443151, 1.53726, -2.40399]
        q_ij[:, 4, 0] = [0.00852515, 1.88625, -2.64651]
        q_ij[:, 5, 0] = [-0.123100, 0.983854, -1.29188]
        q_ij[:, 6, 0] = [-0.0376917, 0.187708, -0.235359]
        #   ---------- Columns ----------
        q_ij[:, 0, 1] = [-0.00189096, 0.000637430, 0.00157383]
        q_ij[:, 1, 1] = [0.00981029, 0.0409220, 0.00908004]
        q_ij[:, 2, 1] = [0.732647, 0.0539796, -0.665773]
        q_ij[:, 3, 1] = [-1.59927, -0.500870, 1.86375]
        q_ij[:, 4, 1] = [1.54047, 0.692547, -2.05390]
        q_ij[:, 5, 1] = [-0.707187, -0.374173, 1.01287]
        q_ij[:, 6, 1] = [0.125276, 0.0721572, -0.186466]

        # refractive index correction of asymmetry parameter
        c_g = np.zeros(4).reshape(2, 2)
        c_g[:, 0] = [0.96025050, 0.42918060]
        c_g[:, 1] = [0.94179149, -0.21600979]
        # correction for absorption
        s = [1.00014, 0.666094, -0.535922, -11.7454, 72.3600, -109.940]
        u = [-0.213038, 0.204016]

        # -------- selector for plates or columns
        if ar > 1.:
            col_pla = 1  # columns
        else:
            col_pla = 0  # plates & compacts

        # ------------ Size parameters -------------------

        # --- absorption size parameter (Fig. 4, box 1)
        Chi_abs = mi / wl * V / Area
        # ----- scattering size parameter (Fig. 7, box 1)
        Chi_scat = 2. * np.pi * np.sqrt(Area / np.pi) / wl

        # ------------ SINGLE SCATTERING ALBEDO ----------

        if Chi_abs > 0:
            w_1 = 1. - a[0] * (1. - np.exp(-Chi_abs * a[1]))  # for AR=1 (Fig. 4, box 2)
            l = np.zeros(nc1)
            for i in range(nc2): l[:] += c_ij[:, i, col_pla] * np.log10(ar) ** i  # (Fig. 4, box 3)
            D_w = l[0] * np.exp(-(np.log(Chi_abs) - l[2]) ** 2 / (2. * l[1] ** 2)) / (
                        Chi_abs * l[1] * np.sqrt(2. * np.pi))  # (Fig. 4, box 3)
            w = w_1 + D_w  # (Fig. 4, box 4)
        else:
            w = 1.

        # --------------- ASYMMETRY PARAMETER ------------

        # diffraction g
        g_diffr = b_gdiffr[0] * np.exp(b_gdiffr[1] * np.log(Chi_scat)) + b_gdiffr[2]  # (Fig. 7, box 2)
        g_diffr = max([g_diffr, 0.5])

        # raytracing g at 862 nm
        g_1 = 0.
        for i in range(len(p_a_eq_1)): g_1 += p_a_eq_1[i] * delta ** i  # (Fig. 7, box 3)

        p_delta = np.zeros(nq1)
        for i in range(nq2): p_delta += q_ij[:, i, col_pla] * np.log10(ar) ** i  # (Fig. 7, box 4)
        Dg = 0.
        for i in range(nq1): Dg += p_delta[i] * delta ** i  # (Fig. 7, box 4)
        g_rt = 2. * (g_1 + Dg) - 1.  # (Fig. 7, box 5)

        # --------- refractive index correction of asymmetry parameter (Fig. 7, box 6)
        epsilon = c_g[0, col_pla] + c_g[1, col_pla] * np.log10(ar)
        mr1 = 1.3038  # reference value @ 862 nm band
        C_m = abs((mr1 - epsilon) / (mr1 + epsilon) * (mr + epsilon) / (
                    mr - epsilon))  # abs function added according to corrigendum to the original paper

        # ---- correction for absorption (Fig. 7, box 7)
        if Chi_abs > 0:
            C_w0 = 0.
            for i in range(len(s)): C_w0 += s[i] * (1. - w) ** i
            k = np.log10(ar) * u[col_pla]
            C_w1 = k * w - k + 1.
            C_w = C_w0 * C_w1
        else:
            C_w = 1.

        # raytracing g at required wavelength
        g_rt_corr = g_rt * C_m * C_w  # (Fig. 7, box 9)

        # ------ Calculate total asymmetry parameter and check g_tot <= 1 (Fig. 7, box 9)
        g_tot = 1. / (2. * w) * ((2. * w - 1.) * g_rt_corr + g_diffr)
        g_tot = min([g_tot, 1.])

        absXS = ((1 - (np.exp(-4 * np.pi * mi * V / Area * wavelengths[i]))) * Area)
        X = (2 * np.pi * Reff) / wl

        # append values to output lists
        X_list.append(X)
        Chi_abs_list.append(Chi_abs)
        SSA_list.append(w)
        Assy_list.append(g_tot)
        absXS_list.append(absXS)
 
    if plots:

        plt.figure(1)
        plt.plot(wavelengths, SSA_list, label='{}x{}'.format(r, depth)), plt.xlim(0.3,1.4), plt.ylabel('SSA'), plt.xlabel(
            'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)

        if savefigs:
            plt.savefig(str(savepath+'SSA_{}x{}.jpg'.format(r,depth)))

        plt.figure(2)
        plt.plot(wavelengths, Assy_list, label='{}x{}'.format(r, depth)), plt.xlim(0.3,1.4), plt.ylabel('Assymetry Parameter'), plt.xlabel(
            'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
        
        if savefigs:
            plt.savefig(str(savepath+'AssymetryParam_{}x{}.jpg'.format(r,depth)))


        plt.figure(3)
        plt.plot(wavelengths, absXS_list, label='{}x{}'.format(r, depth)), plt.xlim(0.3,1.4), plt.ylabel(
            'Absorption Cross Section'), plt.xlabel('Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
        
        if savefigs:
            plt.savefig(str(savepath+'AbsXS_{}x{}.jpg'.format(r,depth)))

        plt.figure(4)
        plt.plot(wavelengths, X_list, label='{}x{}'.format(r, depth)), plt.ylabel('Size Parameter X'), plt.xlabel(
            'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
     
        if savefigs:
            plt.savefig(str(savepath+'X_SizeParam_{}x{}.jpg'.format(r,depth)))

    if report_dims:
        print('cell diameter = ', np.round(diameter, 2), ' (micron)')
        print('cell length = ', depth, ' (micron)')
        print('aspect ratio = ', ar)
        print('cell volume = ', np.round(V, 2), ' (micron^3)')
        print('Effective radius = ', np.round(Reff, 2), ' (micron)')
        print("projected area = ", np.round(Area, 2))
        print()  # line break


    return Assy_list, SSA_list, absXS_list, depth, r, Chi_abs_list, Reff, X_list



def calc_optical_params_MIE(basePath, cell_radius, cell_density, reals, imags, wavelengths, plots=True, savefigs = True):
    """
    cell radius is in microns
    density is in kg/m3
    reals is generated by preprocess_RI()
    imags is generated by preprocess_RI()
    wavelengths is generated by preprocess_RI()
    name is the component of the filename that identifies the specific model run
    plots is a Boolean that toggles plotting of the oprical params

    """

    qqabs = np.zeros(len(wavelengths))
    qqsca = np.zeros(len(wavelengths))
    qqext = np.zeros(len(wavelengths))
    MAC = np.zeros(len(wavelengths))
    assym = np.zeros(len(wavelengths))

    for i in np.arange(0,len(wavelengths),1):
    # Use K to drive mie solver
        
        m = reals[i][0]
        k = imags[i][0]

        X = 2*np.pi*cell_radius/wavelengths[i]

        qext, qsca, qback, g = mie(complex(str("{}-{}j".format(m,k))),X)

        abs_coeff = 4 * np.pi * (imags[i] / wavelengths[i]*1000)

        qabs = qext - qsca
        
        qqabs[i]=qabs*np.pi*cell_radius**2 # calculate cross section from efficiency
        qqsca[i]=qsca*np.pi*cell_radius**2 # calculate cross section from efficiency
        qqext[i]=qext*np.pi*cell_radius**2 # calculate cross section from efficiency

        MAC[i] = abs_coeff/cell_density # calculate mass absorption coefficient from cross section
        assym[i] = g


    ss_alb = qqsca/qqext


    if plots:

        plt.figure(1)
        plt.plot(wavelengths,qqabs,'b',label='absorption cross section')
        plt.plot(wavelengths,qqsca,'g',label='scattering cross section')
        plt.plot(wavelengths,qqext,'r', label= 'extinction cross section')
        plt.xlim(0.3,2.5)
        
        plt.ylabel(r"Cross Section ($\mu$m$^2$)")
        plt.xlabel('Wavelength (nm)')
        plt.legend(loc='best')

        if savefigs:
            plt.savefig(str(basePath + "crosssections{}.jpg".format(name)),dpi=150)
            plt.show()
        else:
            plt.show()

        return assym, ss_alb



def net_cdf_updater(savepath, Assy_list, SSA_list, MAC, depth, r, density):
    
    """

    Function collates optical properties into NetCDF format that can be accessed by RT model

    params: 
    filepath: destination path for NetCDFs
    Assy list: list of assymetry parameter per unit wavelength
    SSA_list: list of siungle scattering albedo per unit wavelength
    MAC: list of mass absorption coefficient per unit wavelength
    depth: cell length (microns)
    r: cell radius (microns)
    density: density of cell (default 1400 kg m-3)

    returns: None

    """

    algfile = pd.DataFrame()
    algfile['asm_prm'] = np.squeeze(Assy_list)
    algfile['ss_alb'] = np.squeeze(SSA_list)
    algfile['ext_cff_mss'] = MAC
    algfile = algfile.to_xarray()
    algfile.attrs['medium_type'] = 'air'
    algfile.attrs['description'] = 'Optical properties for algal cell: cylinder of radius {}um and length {}um'.format(
        str(r), str(depth))
    algfile.attrs['psd'] = 'monodisperse'
    algfile.attrs['radius_um'] = r
    algfile.attrs['side_length_um'] = depth
    algfile.attrs['density_kg_m3'] = density
    algfile.attrs[
        'origin'] = 'Optical properties derived from geometrical optics calculations (algae_go.py) with empirically derived MAC'
    algfile.to_netcdf(str(savepath + 'RealPhenol_algae_geom_{}_{}.nc'.format(str(r), str(depth))), mode='w')

    return