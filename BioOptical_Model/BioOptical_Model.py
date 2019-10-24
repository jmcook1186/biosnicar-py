#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 15:24:05 2018

@author: joe

This code predicts the refractive index and mass absorption coefficient of cylindrical algal cells whose pigment
profile is known. There are several ways to run this script, offering blends of theoretical and empirical
pigment data. The user must first select whether the mass absorption coefficient of the algal cel is to be
calculated from the absorption coefficients of the individual pigments and user defined pigment profile, or whether
the mass absorption coefficient should be loaded in from an external file. The latter option is used to provide
values derived from empirical measurements. This is achieved by setting either calc_MAC or load_MAC to true and the other
to False in the function call.

The pigment profile can then be provided as actual mass of pigment (ng) per cell along with the mass of the cell (ng)
or provided as a mass fraction (% total cellular dry weight). This is achieved by setting pig_mass or pig_frac to True
and the other to False in the function call. There is also an option to apply a correction for packaging effects in
the cell. If the MAC is loaded in from an external file, setting apply_packaging_correction = True loads in a different
csv file where the correction has already been applied. If the MAC is being calculated on-the-fly then the correction
factor per wavelength is loaded in for an element-wise multiplication of the MAC and the correction factor.

The script will then return the cell MAC (m2/kg) and refractive index (dimensionless). The refractive index is
predicted using a mixing model adapted from Pottier et al (2005), Dauchet et al (2015) and Cook et al (2017).

Both the MAC and refractive index are required for the script "Algae_GO.py" to calculate single scattering optical properties
that enable the algae to be incorporated as an impurity in BioSNICAR_GO. The paths in this script are set so that files
are automatically saved in directories accessible to Algae_Go.py. Any changes to paths should also be updated in
Algae_GO.py.

Together, this script and Algae_GO.py populate BioSNICAR_GO's lookup library of netCDF algal optical properties.

*****************************************************************************************************************
*** NB STRONGLY RECOMMEND EMPIRICAL MAC IS LOADED IN FROM FILE AND CORRECTED FOR PACKAGING EFFECTS AS DEFAULT ***
********************* ALSO SUGGEST PROVIDING MEASURED MASS OF PIGMENT IN MG WHERE POSSIBLE **********************
*****************************************************************************************************************

INPUT PARAMETERS:

load_MAC: if True, MAC loaded from external file
apply_packaging_correction: correct the MAC for packaging effects using a wavelength dependent correction factor (CW)
calc_MAC: If true, MAC calculated theoretically
calc_k: if True, script calculates refractive index
pig_mass: if True, pigment data should be provided as absolute mass of pigment per cell (mg)
pig_frac: if True, pigment data should be provided as mass fraction of pigment per cell (% total dry weight)
Pottier: If true, use Pottier et al's (2005) equation to predict the imaginary refractive index
Cook: If true, use Cook et al's (2019) updated equation to predict the imaginary refractive index
cell_dm_weight: provide dry weight of cell in ng
chla,chlb,ppro,psyn,purp: either mass (mg) or mass fraction of each pigment in cell
Xw: water fraction, set to 0.8 as default
density: density of dry cellular material in kg/m3, set to 1400 as default
nm: real part of refractive index, set to 1.4 as default
savefiles: if True, the data will be saved as a csv file to path specified in savefilename
savefilename: path to save files
plot_title: title for plots
plot_optical_props: if true, k and MAC plotted in ipython console
plot_pigment_MACs: if true, the mass absorption coefficient of each pigment is plotted to ipython console

OUTPUTS:

k_list: imaginary refractive index per unit wavelength from 300-5000 nm in at 10nm resolution
MAC: mass absorption coefficient per unit wavelength from 300-5000 nm at 10nm resolution
data: pandas dataframe containing nm, MAC, k per unit wavelength from 300 - 5000 nm at 10nm resolution

"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

def bio_optical(load_MAC = True, apply_packaging_correction=True, calc_MAC = False, calc_k = True, pig_mass = True,
                pig_frac = False, Pottier = False, Cook = True, cell_dm_weight = 0.82, chla = 0.01, chlb = 0.00066,
                ppro = 0.01, psyn = 0, purp = 0.068, Xw = 0.8, density= 1400, nm = 1.4, savefiles = False,
                savepath = "path", savefilename = "name", plot_optical_props = True, plot_pigment_MACs = True):

    data = pd.DataFrame() # set up dataframe for storing optical properties

    if load_MAC: # choose this option to load an empirically derived MAC from file

        if apply_packaging_correction: # optionally apply correction for packaging effect
            MAC = pd.read_csv('/home/joe/Code/BioSNICAR_GO/phenol_mac_packaging_corr.csv',header=None,names=['MAC'])
            MAC = MAC[0:4695] # subsample to appropriate resolution for snicar
            MAC = MAC[0:-1:10]
            data['MAC'] = MAC['MAC'].dropna() # drop NaNs and save to dataframe

        else:
            MAC = pd.read_csv('/home/joe/Code/BioSNICAR_GO/Empirical_MAC.csv',header=None,names=['MAC'])
            MAC = MAC[0:4695] # subsample to appropriate resolution for snicar
            MAC = MAC[0:-1:10]
            data['MAC'] = MAC['MAC'].dropna() # drop NaNs and save to dataframe

    if calc_MAC or calc_k: # choose this option to calculate MAC or k theoretically
        # set up empty lists
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
    
        # Read in wavelength dependent in vivo absorption coefficients for each pigment
        # and define wavelength range
        
        WL = np.arange(300,5000,10)
        WLmeters = [i*1e-9 for i in (WL)]
    
        with open('/home/joe/Code/chlorophyll-a.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea1.append(cellf)
    
        with open('/home/joe/Code/chlorophyll-b.csv') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea2.append(cellf)
    
        with open('/home/joe/Code/Photoprotective_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea3.append(cellf)
    
        with open('/home/joe/Code/Photosynthetic_carotenoids.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)*1e-6
                    Ea4.append(cellf)
    
        with open('/home/joe/Code/BioSNICAR_GO/phenol_mac_correction.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    Ea5.append(cellf)
    
        with open('/home/joe/Code/water_RI.csv')as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    WatRI.append(cellf)
        # read in correction factor for packaging effect (from Chris W paper)
        with open('/home/joe/Code/BioSNICAR_GO/phenol_mac_correction.csv')as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                for cell in row:
                    cellf = float(cell)
                    MACcorr.append(cellf)

        MACcorr = MACcorr[19:-1] # remove wavelengths <300nm
        for i in np.arange(len(MACcorr),750,1): # ensure list length = 750
            MACcorr.append(1)

       # extend spectral data down to 300nm by padding with the value at 350nm and
       # small but non-zero value above 750 nm (avoid divide by zero errors)    
        Ea1 = [Ea1[0] for _ in range(50)] + Ea1
        Ea2 = [Ea2[0] for _ in range(50)] + Ea2
        Ea3 = [Ea3[0] for _ in range(50)] + Ea3
        Ea4 = [Ea4[0] for _ in range(50)] + Ea4
        
        # extend data with zeros at nonabsoring wavelengths to 5000 nm
        for i in np.arange(751,5000,1):
            Ea1.append(0)
            Ea2.append(0)
            Ea3.append(0)
            Ea4.append(0)
            Ea5.append(0)
            MACcorr.append(1)
    
        # downsample to match SNICAR resolution) 
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

        # copy water RI over zeros at non-absorbing wavelengths
        Ea1n[44:-1] = WatRIn[44:-1]
        Ea2n[44:-1] = WatRIn[44:-1]
        Ea3n[44:-1] = WatRIn[44:-1]
        Ea4n[44:-1] = WatRIn[44:-1]
        Ea5n[44:-1] = WatRIn[44:-1]
            
    if calc_MAC: # if MAC is calculated on-the-fly using theoretical approach
        if pig_frac:
        # if data was provided as mass fraction convert percentage dw pigment to actual mass of pigment per cell (mg)
            chla_w = chla * cell_dm_weight*1e-6
            chlb_w = chlb * cell_dm_weight*1e-6
            ppro_w = ppro * cell_dm_weight*1e-6
            psyn_w = psyn * cell_dm_weight*1e-6
            purp_w = purp * cell_dm_weight*1e-6

        elif pig_mass:
        # If data was provided in units of mg pigment/cell, read in from file.
            chla_w = chla
            chlb_w = chlb
            ppro_w = ppro
            psyn_w = psyn
            purp_w = purp

        # Multiply mass of each pigment (mg) by absorption coefficient (m2/mg)
        # at each wavelenth giving units of m2 per cell
        EW1m = [a * chla_w for a in Ea1n]
        EW2m = [a * chlb_w for a in Ea2n]
        EW3m = [a * ppro_w for a in Ea3n]
        EW4m = [a * psyn_w for a in Ea4n]
        EW5m = [a * purp_w for a in Ea5n]

        EWWm = [sum(x) for x in zip(EW1m,EW2m,EW3m,EW4m,EW5m)] # Sum all pigments (m2)
        MAC = [(i/cell_dm_weight)*1e12 for i in EWWm] #normalise from cell to mass (x 1e12 to get m3/kg)
        data['MAC'] = MAC # save to dataframe

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
        for i in np.arange(0,len(WL),1):
            if Pottier: # Pottier et al (2005) equation
                k = (((1 - Xw) / Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) #original Pottier equation
            elif Cook: # Cook et al (2019) update
                k = (Xw * WatRIn[i]) + ((1 - Xw) * (WLmeters[i]/np.pi*4) * density * EWW[i]) # Cook (2018) updated version
            k_list.append(k)
            real_list.append(nm)

        data['Imag'] = k_list # append imaginary RI to dataframe
        data['Real'] = real_list # append real RI to dataframe

    if savefiles: # optional save dataframe to csv files
        data.to_csv(str(savepath+'{}_Dataset.csv'.format(savefilename)))
        data['Imag'].to_csv(str(savepath+'{}_KK.csv'.format(savefilename)),header=None,index=False)
        data['MAC'].to_csv(str(savepath+'{}_MAC.csv'.format(savefilename)),header=None,index=False)
        data['Real'].to_csv(str(savepath+'{}_Real.csv'.format(savefilename)),header=None,index=False)

    if plot_optical_props:
        plt.figure(figsize=(10,15))
        plt.subplot(2,1,1)
        plt.plot(WL[0:220],MAC[0:220]),plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('MAC (kg/m^3)',fontsize=16)
        plt.tight_layout()

        plt.subplot(2,1,2)
        plt.plot(WL[0:220],k_list[0:220]), plt.xlim(300,750)
        plt.xticks(fontsize=16), plt.yticks(fontsize=16)
        plt.xlabel('Wavelength',fontsize=16),plt.ylabel('K (dimensionless)',fontsize=16)
        plt.tight_layout()

    if plot_pigment_MACs:
        plt.figure(figsize=(12, 12))
        plt.plot(WL, Ea1n, label='Chlorophyll a')
        plt.plot(WL, Ea2n, label='Chlorophyll b')
        plt.plot(WL, Ea3n, label='Secpndary carotenoids')
        plt.plot(WL, Ea4n, label='Primary carotenoids')
        plt.plot(WL, Ea5n, label='Purpurogallin')
        plt.xlabel('Wavelengths nm')
        plt.ylabel('In vivo mass absorption coefficient (m2/kg)')
        plt.xlim(300, 750)

        plt.legend(loc='best')
        plt.show()

        return k_list, real_list, MAC, data


# NB pigment data is provided here in units of mg per cell
k_list, real_list, MAC, data = bio_optical(
        load_MAC= False,
        apply_packaging_correction=True,
        calc_MAC = True,
        calc_k = True,
        pig_mass = False,
        pig_frac = True,
        Pottier = False,
        Cook = True,
        cell_dm_weight= 1.89,
        chla = 0.028,
        chlb = 0.013,
        ppro = 0,
        psyn = 0.063,
        purp = 0,
        Xw = 0.8, 
        density= 1400, 
        nm = 1.4, 
        savefiles = False,
        savepath = '/home/joe/Code/BioSNICAR_GO/',
        savefilename = "snw_alg1_20_7_2019",
        plot_optical_props = True,
        plot_pigment_MACs = True)

