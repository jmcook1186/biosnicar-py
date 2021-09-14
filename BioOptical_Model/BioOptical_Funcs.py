#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jmcook1186, lcvl41

Contains functions relating to bio-optical model components of BioSNICAR_GO_py
Called from BioOptical_driver.py

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.signal import savgol_filter
from miepython import mie
from plotnine import ggplot, aes, geom_line
plt.style.use('seaborn')


def bioptical_calculations(ACS_calculated, ACS_file, biovolume, density, xw, 
                           cell_volume, wvl, packaging_correction, 
                           pigment_dir, pigments_data, n_algae, k_water,
                           smooth, window_size, poly_order, smoothStart, 
                           smoothStop, plot_optical_props, 
                           plot_pigment_ACSs, savefiles, savefilename,
                           saveplots, savepath):
    
    ##################
    ## Initialization
    ################## 
    data = pd.DataFrame() # storing n,k,ACS
    abs_cff_pigments = pd.DataFrame(index = wvl*1000) # storing pigment MACs
    wvl = wvl # µm
    
    ##################
    ## ACS calculation
    ##################      
    # open mass absorption coefficients (m2/mg) for each algal pigment
    # from a dictionary ('pigments_data')
    
    if ACS_calculated: 
        abs_coeff = 0 
        for key,value in pigments_data.items(): # key is pigment name, value is abs coeff in m2/mg
            abs_pigm = (np.array(pd.read_csv(key))).flatten() # store abs coeff in m2/mg 
            abs_cff_pigments[str(key.split(pigment_dir,1)[1])[0:-4]] = abs_pigm # store pigment MACs for plotting
            conc = value # store cellular pigment concentration in mg/µm3 or ng/cell
            abs_coeff = abs_coeff + conc * abs_pigm # in m2/µm3 cells or 10e6 m2/cell
        if biovolume:
            ACS = abs_coeff/density # from m2/µm3 to m2/kg wet biomass, dividing by density in kg/µm3
        else:
            ACS = abs_coeff/1000000 # from 10e6 m2/cell to m2/cell
        
        # apply packaging correction
        if packaging_correction:
            print('not ready yet')
    
    else:
        ACS = (np.array(pd.read_csv(ACS_file))).flatten()
    
    ##################
    ## k calculation
    ##################
    
    if biovolume:
        k = ACS * density * wvl / 4 / np.pi * 10**(12) # units: abs coeff (converted from m2/µm3 to µm2/µm3 so *10**12) * wvl (µm)
        n=n_algae
    else:         
        k = (xw * k_water) + ((1 - xw) * (wvl/(np.pi*4)) / cell_volume * ACS * 10**(12)) # units: ACS (converted from m2/cell to µm2/cell) / volume cell (µm3) * wvl (µm)
        n=n_algae
    
    ###############################
    ## ACS, k storage and rescaling
    ###############################
    
    wvl_rescaled_BioSNICAR = wvl[::10]
    ACS_rescaled_BioSNICAR = ACS[::10]
    k_rescaled_BioSNICAR = k[::10]
    n_rescaled_BioSNICAR = n[::10]

    data['wvl'] = wvl_rescaled_BioSNICAR
    data['ACS'] = ACS_rescaled_BioSNICAR # in m2/kg dry or wet biomass
    data['k'] = k_rescaled_BioSNICAR # unitless
    data['n'] = n_rescaled_BioSNICAR # unitless
    

    ################################
    ## optional: plotting and saving
    ################################

    
    # apply optional smoothing filter with user defined window width and polynomial order
    if smooth:
        yhat = savgol_filter(ACS[smoothStart:smoothStop], window_size, poly_order) # window size 51, polynomial order 3
        ACS[smoothStart:smoothStop] = yhat
    
    # optionally save files to savepath
    if savefiles: # optional save dataframe to csv files
        data['k'].to_csv(str(savepath+'{}_k.csv'.format(savefilename)),header=None,index=False)
        data['ACS'].to_csv(str(savepath+'{}_ACS.csv'.format(savefilename)),header=None,index=False)
        data['n'].to_csv(str(savepath+'{}_n.csv'.format(savefilename)),header=None,index=False)

    # optionally plot figures to interative window
    if plot_optical_props:
        plt.figure(figsize=(10,15))
        plt.subplot(2,1,1)
        plt.plot(wvl[100:700]*1000,ACS[100:700])
        plt.xticks(fontsize=24), plt.yticks(fontsize=24)
        plt.xlabel('Wavelength (nm)',fontsize=24),plt.ylabel(str('ACS (m$2$ cell$^{-1}$ or m$^2$ µm$3$))'),fontsize=24)
        plt.tight_layout()
    
        plt.subplot(2,1,2)
        plt.plot(wvl[100:700]*1000,k[100:700])
        plt.xticks(fontsize=24), plt.yticks(fontsize=24)
        plt.xlabel('Wavelength (nm)',fontsize=24),plt.ylabel('k',fontsize=24)
        plt.tight_layout()
        
        # optionally save plots to savepath
        if saveplots:
            plt.show()
            plt.savefig(str(savepath+"/ACSandK.png"))
        else:
            plt.show()
        
    if plot_pigment_ACSs:
        ax=abs_cff_pigments.replace(0, np.nan).plot(use_index=True,xlim=(300, 800))
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Mass absorption coefficient (m$^2$ mg$^{-1}$)")
        abs_cff_pigments.plot()
    
        #if saveplots:
            #plt.show()
            #plt.savefig(str(savepath+"/PigmentACSs.png"))
        
        #else:
            #plt.show()
    

    return wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, n_rescaled_BioSNICAR, ACS, ACS_rescaled_BioSNICAR, data, abs_cff_pigments


#%%

def ssp_calculations(GO, Mie, savepath, r, L, wvl, n_algae, k_algae, plots, savefigs, figname, report_dims):
    
    if GO:     
        # set up lists
        SSA_list = []
        Assy_list = []
        absXS_list = []
        Chi_abs_list = []
        X_list = []
        
        # calculate cylinder dimensions
        diameter = 2 * r
        V = L * (np.pi * r ** 2)  # volume of cylinder in µm3
        Reff = (V / ((
                                 4 / 3) * np.pi)) ** 1 / 3  # effective radius (i.e. radius of sphere with equal volume to real cylinder)
        Area_total = 2 * (np.pi * r ** 2) + (2 * np.pi * r) * (
            L)  # total surface area - 2 x ends plus circumference * length
        Area = Area_total / 4  # projected area
        ar = diameter / L  # aspect ratio
        delta = 0.3
        
        # start geometric optics calculation per wavelength
        for i in np.arange(0, len(wvl), 1):
            mr = n_algae[i]  # load real RI per wavelength
            mi = k_algae[i]  # load imaginary RI per wavelength
            wl = wvl[i]  # loadd wavelength
        
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
        
            absXS = ((1 - (np.exp(-4 * np.pi * mi * V / Area * wvl[i]))) * Area)
            X = (2 * np.pi * Reff) / wl
        
            # append values to output lists
            X_list.append(X)
            Chi_abs_list.append(Chi_abs)
            SSA_list.append(w)
            Assy_list.append(g_tot)
            absXS_list.append(absXS)
            # convert to np arrays to return by the function
            assym=np.array(Assy_list).flatten()
            ss_alb=np.array(SSA_list).flatten()
            
        
         
        if plots:
        
            plt.figure(1)
            plt.plot(wvl, SSA_list, label='{}x{}'.format(r, L)), plt.xlim(0.3,1.4), plt.ylabel('SSA'), plt.xlabel(
                'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
        
            if savefigs:
                plt.savefig(str(savepath+'SSA_{}x{}.jpg'.format(r,L)))
        
            plt.figure(2)
            plt.plot(wvl, Assy_list, label='{}x{}'.format(r, L)), plt.xlim(0.3,1.4), plt.ylabel('Assymetry Parameter'), plt.xlabel(
                'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
            
            if savefigs:
                plt.savefig(str(savepath+'AssymetryParam_{}x{}.jpg'.format(r,L)))
        
        
            plt.figure(3)
            plt.plot(wvl, absXS_list, label='{}x{}'.format(r, L)), plt.xlim(0.3,1.4), plt.ylabel(
                'Absorption Cross Section'), plt.xlabel('Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
            
            if savefigs:
                plt.savefig(str(savepath+'AbsXS_{}x{}.jpg'.format(r,L)))
        
            plt.figure(4)
            plt.plot(wvl, X_list, label='{}x{}'.format(r, L)), plt.ylabel('Size Parameter X'), plt.xlabel(
                'Wavelength (um)'), plt.grid(False), plt.legend(loc='best', ncol=2)
         
            if savefigs:
                plt.savefig(str(savepath+'X_SizeParam_{}x{}.jpg'.format(r,L)))
        
        if report_dims:
            print('cell diameter = ', np.round(diameter, 2), ' (micron)')
            print('cell length = ', L, ' (micron)')
            print('aspect ratio = ', ar)
            print('cell volume = ', np.round(V, 2), ' (micron^3)')
            print('Effective radius = ', np.round(Reff, 2), ' (micron)')
            print("projected area = ", np.round(Area, 2))
            print()  # line break
        
        
    if Mie: 

        X = 2*np.pi*r/wvl # unitless
        qext, qsca, qback, g = mie(n_algae-1j*k_algae,X)        
        qabs = qext - qsca
        qqabs=qabs*np.pi*r**2 # calculate cross section from efficiency
        qqsca=qsca*np.pi*r**2 # calculate cross section from efficiency
        qqext=qext*np.pi*r**2 # calculate cross section from efficiency
        assym = g
        ss_alb = qqsca/qqext

        if plots:
    
            plt.figure(1)
            plt.plot(wvl,qqabs,'b',label='absorption cross section')
            plt.plot(wvl,qqsca,'g',label='scattering cross section')
            plt.plot(wvl,qqext,'r', label= 'extinction cross section')
            plt.xlim(0.2,2.5)
            
            plt.ylabel(r"Cross Section ($\mu$m$^2$)")
            plt.xlabel('Wavelength (nm)')
            plt.legend(loc='best')
    
            if savefigs:
                plt.savefig(str(savepath + figname + ".jpg"),dpi=150)
                plt.show()
            else:
                plt.show()
        
    return assym, ss_alb
    


def net_cdf_updater(GO, Mie, savepath, filename, wvl, g, ssa, ACS, L, r, density, information):

    algfile = pd.DataFrame()
    algfile['asm_prm'] = np.squeeze(g)
    algfile['ss_alb'] = np.squeeze(ssa)
    algfile['ext_cff_mss'] = ACS
    algfile = algfile.to_xarray()
    algfile.attrs['medium_type'] = 'air'
    if GO: 
        algfile.attrs['description'] = 'Optical properties for glacier algal cell: cylinder of radius {}um and length {}um'.format(
        str(r), str(L))
        algfile.attrs['side_length_um'] = L
    if Mie: 
        algfile.attrs['description'] = 'Optical properties for snow algal cell: sphere of radius {}um'.format(
        str(r))
    algfile.attrs['psd'] = 'monodisperse'
    algfile.attrs['radius_um'] = r
    algfile.attrs['density_kg_µm3'] = density
    algfile.attrs['wvl'] = wvl
    algfile.attrs['information'] = information
    algfile.to_netcdf(str(savepath + filename + '.nc'), mode='w')

    return