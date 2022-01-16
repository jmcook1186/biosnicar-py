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


def bioptical_calculations(ACS_calculated, ACS_file, ACS_loaded_invivo,
                           ACS_loaded_reconstructed, biovolume, biomass, 
                           cellular, density_wet, density_dry,
                           dir_base, cell_vol, wvl, 
                           packaging_correction_SA,
                           packaging_correction_GA, 
                           pigment_dir, pigments_data, n_algae, k_water,
                           smooth, window_size, poly_order, smoothStart, 
                           smoothStop, plot_optical_props, 
                           plot_pigment_ACSs, savefiles, savefilename,
                           saveplots, savepath):
    
    #################
    ## Initialization
    ################# 

    data = pd.DataFrame() # storing n,k,ACS
    abs_cff_pigments = pd.DataFrame(index = wvl * 1000) # storing pigment MACs
    wvl = wvl # µm
    xw = 0.59 * density_wet / 1000 # conversion mass to volume water fraction
    
    ##################
    ## ACS calculation
    ##################     

    if ACS_calculated: 
    # open mass absorption coefficients (m2/mg) for each algal pigment
    # from a dictionary.key is pigment name, value is abs coeff in m2/mg
        abs_coeff = 0 
        for key,value in pigments_data.items(): 
            abs_pigm = (np.array(pd.read_csv(key))).flatten() # m2/mg 
            abs_cff_pigments[str(key.split(pigment_dir,1)[1])[0:-4]]=abs_pigm 
            conc = value # intracellular conc in ng/µm3, ng/cell, or ng/ng
            abs_coeff = abs_coeff + conc*abs_pigm # 10e6 m2/µm3,m2/cell,m2/ng
        
        if biovolume: 
            ACS = abs_coeff/1000000 # correct the 10e6 
        if cellular: 
            ACS = abs_coeff/1000000 # correct the 10e6 
        if biomass: 
            ACS = abs_coeff/1000 # m2/mg for biomass option
    
    elif ACS_loaded_reconstructed:
        ACS = (np.array(pd.read_csv(ACS_file))).flatten() #m2/mg, um3 or cell
        if packaging_correction_SA: # !! only from 300nm, rest set to 0
            pckg_SA = np.loadtxt(dir_base + 'Data/pigments/pckg_SA.csv')
            abs_coeff = abs_coeff * pckg_SA[0:-1] 
        if packaging_correction_GA: # !! only from 300nm, rest set to 0
            pckg_GA = np.loadtxt(dir_base + 'Data/pigments/pckg_GA.csv')
            abs_coeff = abs_coeff * pckg_GA[0:-1]
            
    elif ACS_loaded_invivo:
        ACS = (np.array(pd.read_csv(ACS_file))).flatten() #m2/mg, um3 or cell
        
    ################
    ## k calculation
    ################ 
    
    k_water_alg = k_water
    k_water_alg[0:600] = 0 
    
    if cellular: # ACS in m2 / cell
    # units: ACS (m2/cell to µm2/cell) / cell volume (um3/cell) * wvl (µm)
        k = xw * k_water_alg + ACS * 10**(12) * wvl / (np.pi * 4) / cell_vol
        n = n_algae 
        
    if biovolume: # ACS in m2 / µm3
    # units: ACS (m2/µm3 to µm2/µm3) * wvl (µm)
        k = xw * k_water_alg + ACS * 10**(12) * wvl / (np.pi * 4)  
        n = n_algae
        
    if biomass: # ACS in m2 / dry mg
    # units: ACS (m2/mg to m2/kg) * density (kg m3) * wvl (µm to m)
        k = xw * k_water + ACS * density_dry * wvl / (np.pi * 4)
        n = n_algae
    
    ###############################
    ## ACS, k storage and rescaling
    ###############################
    
    # rescaling variables to BioSNICAR resolution (10nm)
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

    
    # apply optional smoothing filter with user defined window width 
    # and polynomial order
    if smooth:
        yhat = savgol_filter(ACS[smoothStart:smoothStop], window_size, 
                             poly_order) # window size 51, polynomial order 3
        ACS[smoothStart:smoothStop] = yhat
    
    # optionally save files to savepath
    if savefiles: # optional save dataframe to csv files
        data['k'].to_csv(str(savepath+'{}_k.csv'.format(savefilename)),
                         header=None,index=False)
        data['ACS'].to_csv(str(savepath+'{}_ACS.csv'.format(savefilename)),
                           header=None,index=False)
        data['n'].to_csv(str(savepath+'{}_n.csv'.format(savefilename)),
                         header=None,index=False)

    # optionally plot figures to interative window
    if plot_optical_props:
        plt.figure(figsize=(10,15))
        plt.subplot(2,1,1)
        plt.plot(wvl[100:600]*1000,ACS[100:600])
        plt.xticks(fontsize=24), plt.yticks(fontsize=24)
        plt.xlabel('Wavelength (nm)',fontsize=24),
        plt.ylabel(str('ACS (m$2$ cell$^{-1}$, m$^2$ µm$3$ or m$^2$ ng$3$ ))'),
                   fontsize=24)
        plt.tight_layout()
    
        plt.subplot(2,1,2)
        plt.plot(wvl[100:600]*1000,k[100:600])
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
        ax=abs_cff_pigments.replace(0, np.nan).plot(use_index=True,
                                                    xlim=(300, 800))
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Mass absorption coefficient (m$^2$ mg$^{-1}$)")
        abs_cff_pigments.plot()
    
        #if saveplots:
            #plt.show()
            #plt.savefig(str(savepath+"/PigmentACSs.png"))
        
        #else:
            #plt.show()
    

    return wvl, wvl_rescaled_BioSNICAR, k, k_rescaled_BioSNICAR, n, \
            n_rescaled_BioSNICAR, ACS, ACS_rescaled_BioSNICAR, data, \
            abs_cff_pigments


#%%

def ssp_calculations(GO, Mie, savepath, r, L, wvl, n_algae, k_algae, plots, 
                     savefigs, figname, report_dims):
    
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
        Reff = (V / ((4 / 3) * np.pi)) ** 1 / 3  # effective radius 
        # (i.e. radius of sphere with equal volume to real cylinder)
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
        
            # ---- g correction for AR != 1 (Also applied to AR=1 as plate) 
            # (Table 3)
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
                # for AR=1 (Fig. 4, box 2)
                w_1 = 1. - a[0] * (1. - np.exp(-Chi_abs * a[1]))  
                l = np.zeros(nc1)
                for i in range(nc2): 
                    # (Fig. 4, box 3)
                    l[:] += c_ij[:, i, col_pla] * np.log10(ar) ** i  
                D_w = l[0] * np.exp(-(np.log(Chi_abs) - l[2]) ** 2 
                      / (2. * l[1] ** 2)) / (Chi_abs * l[1] * 
                      np.sqrt(2. * np.pi))  # (Fig. 4, box 3)
                w = w_1 + D_w  # (Fig. 4, box 4)
            else:
                w = 1.
        
            # --------------- ASYMMETRY PARAMETER ------------
        
            # diffraction g
            g_diffr = b_gdiffr[0] * np.exp(b_gdiffr[1] * np.log(Chi_scat)) + \
            b_gdiffr[2]  # (Fig. 7, box 2)
            g_diffr = max([g_diffr, 0.5])
        
            # raytracing g at 862 nm
            g_1 = 0.
            # (Fig. 7, box 3)
            for i in range(len(p_a_eq_1)): g_1 += p_a_eq_1[i] * delta ** i  
        
            p_delta = np.zeros(nq1)
            # (Fig. 7, box 4)
            for i in range(nq2): 
                p_delta += q_ij[:, i, col_pla] * np.log10(ar) ** i  
            Dg = 0.
            # (Fig. 7, box 4)
            for i in range(nq1): Dg += p_delta[i] * delta ** i  
            g_rt = 2. * (g_1 + Dg) - 1.  # (Fig. 7, box 5)
        
            # --------- ref idx correction of asym parameter (Fig. 7, box 6)
            epsilon = c_g[0, col_pla] + c_g[1, col_pla] * np.log10(ar)
            mr1 = 1.3038  # reference value @ 862 nm band
            # abs function added according to corrigendum to original paper
            C_m = abs((mr1 - epsilon) / (mr1 + epsilon) * (mr + epsilon) / (
                        mr - epsilon))  
            
        
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
        
            # ------ Calc tot asym parameter + check g_tot <= 1 (Fig. 7, box 9)
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
            plt.plot(wvl, SSA_list, label='{}x{}'.format(r, L)), 
            plt.xlim(0.3,1.4), plt.ylabel('SSA'), 
            plt.xlabel('Wavelength (um)'), plt.grid(False), 
            plt.legend(loc='best', ncol=2)
        
            if savefigs:
                plt.savefig(str(savepath+'SSA_{}x{}.jpg'.format(r,L)))
        
            plt.figure(2)
            plt.plot(wvl, Assy_list, label='{}x{}'.format(r, L)), 
            plt.xlim(0.3,1.4), plt.ylabel('Assymetry Parameter'), 
            plt.xlabel('Wavelength (um)'), plt.grid(False), 
            plt.legend(loc='best', ncol=2)
            
            if savefigs:
                plt.savefig(str(savepath+'AssymetryParam_{}x{}.jpg'.
                                format(r,L)))
        
        
            plt.figure(3)
            plt.plot(wvl, absXS_list, label='{}x{}'.format(r, L)), 
            plt.xlim(0.3,1.4), plt.ylabel('Absorption Cross Section'), 
            plt.xlabel('Wavelength (um)'), plt.grid(False), 
            plt.legend(loc='best', ncol=2)
            
            if savefigs:
                plt.savefig(str(savepath+'AbsXS_{}x{}.jpg'.format(r,L)))
        
            plt.figure(4)
            plt.plot(wvl, X_list, label='{}x{}'.format(r, L)), 
            plt.ylabel('Size Parameter X'), plt.xlabel('Wavelength (um)'), 
            plt.grid(False), plt.legend(loc='best', ncol=2)
         
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
        assym = g
        ss_alb = qsca/qext

        if plots:
            qqabs=qabs*np.pi*r**2 # calculate cross section from efficiency
            qqsca=qsca*np.pi*r**2 # calculate cross section from efficiency
            qqext=qext*np.pi*r**2 # calculate cross section from efficiency
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
    


def net_cdf_updater(GO, Mie, savepath, filename, wvl, g, ssa, ACS, L, r, 
                    density, information):

    algfile = pd.DataFrame()
    algfile['asm_prm'] = np.squeeze(g)
    algfile['ss_alb'] = np.squeeze(ssa)
    algfile['ext_cff_mss'] = ACS
    algfile = algfile.to_xarray()
    algfile.attrs['medium_type'] = 'air'
    if GO: 
        algfile.attrs['description'] = \
        'Optical properties for glacier algal cell: cylinder of radius '\
        '{}um and length {}um'.format(str(r), str(L))
        algfile.attrs['side_length_um'] = L
    if Mie: 
        algfile.attrs['description'] = 'Optical properties for snow algal '\
            'cell: sphere of radius {}um'.format(
        str(r))
    algfile.attrs['psd'] = 'monodisperse'
    algfile.attrs['radius_um'] = r
    algfile.attrs['density_kg_m3'] = density
    algfile.attrs['wvl'] = wvl
    algfile.attrs['information'] = information
    algfile.to_netcdf(str(savepath + filename + '.nc'), mode='w')

    return