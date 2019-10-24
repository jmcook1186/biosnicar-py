#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 08:38:44 2018

@author: Joseph Cook, joe.cook@sheffield.ac.uk

This file calculates the optial properties (single scattering albedo, assymetry parameter,
mass absorption coefficient and extinction, scattering and absorption cross sections) for
ice grains shaped as arbitrarily large hexagonal plates or columns. The optical properties
are then saved into netCDF files in the correct format for loading into the BioSNICAR_GO 
model provided in the BioSNICAR_GO repository.

The main function calc_optical_params() is based upon the equations of Diedenhoven et al (2014)
who provided a python script as supplementary material for their paper:

"A flexible parameterization for shortwave optical properties of ice crystals" by 
Bastiaan van Diedenhoven; Andrew S. Ackerman; Brian Cairns; Ann M. Fridlind
accepted for publication in J. Atmos. Sci. (2013)

The original code can be downloaded from:
https://www.researchgate.net/publication/259821840_ice_OP_parameterization

The optical properties are calculated using a parameterization of geometric optics
calculations (Macke et al., JAS, 1996).

The script is divided into three functions. The first is a preprocessing function
that ensures the wavelengths and real/imaginary parts of the refractive index for ice is
provided in the correct waveband and correct spectral resolution to interface with the
BioSNICAR_GO model. The refractive indices are taken from Warren and Brandt 2008.

There are no user defined inouts for the preprocessing function, it can simply be
run as 

reals, imags, wavelengths = preprocess()

The calc_optical_params() fnction takes several inputs. reals, imags and wavelengths
are output by preprocess() and side_length and depth are user defined. These are the two
parameters that control the dimensions of the ice crystals. Side_length is the length
in microns of one side of the hexagnal face of the crystal, depth is the column length
also in microns describing the z dimension. The code then calculates volume, apotherm,
aspect ratio, area etc inside the function. The optical parameters are returned.
Optional plots and printed values for the optical params are provided by setting
plots to true and the dimensions of the crystals can be reported by setting
report_dims to true in the function call.

The final function, net_cdf_updater() is used to dump the optical parameters and
metadata into a netcdf file and save it into the working directory to be used as
a lookup library for the two-stream radiative transfer model BoSNICAR_GO.

The function calls are provided at the bottom of this script in a loop, where the
user can define the range of side lengths and depths to be looped over.

NOTE: The extinction coefficient in the current implementation is 2 for all size parameters 
as assumed in the conventional geometric optics approximation.

UPDATE March 2019: Now includes option to output a three-band version, with bands relevant for integrating
BioSNICAR into MAR (modele atmospherique regionale). These bands are 300-799 nm, 800-1199nm, 1201-2500nm. The
optical properties are averaged over these wavelength ranges and output as three floats to a separate optical
property library: GO_files/Ice_Optical_Properties_3band. To use this option set ThreeBand to True in the functions
calculate_optical_params() and netcdf_updater()

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd

filepath = '/home/joe/Code/BioSNICAR_GO/GO_files/'
filepath3Band = '/home/joe/Code/BioSNICAR_GO/GO_files/Ice_Optical_Properties_3band/' #save path for 3 band version


def preprocess_RI():
    
    wavelengths = [0.210, 0.250, 0.300,
         0.350, 0.390, 0.400, 0.410, 0.420, 0.430, 0.440, 0.450, 0.460, 0.470, 0.480,
         0.490, 0.500, 0.510, 0.520, 0.530, 0.540, 0.550, 0.560, 0.570, 0.580, 0.590,
         0.600, 0.610, 0.620, 0.630, 0.640, 0.650, 0.660, 0.670, 0.680, 0.690, 0.700,
         0.710, 0.720, 0.730, 0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800, 0.810,
         0.820, 0.830, 0.840, 0.850, 0.860, 0.870, 0.880, 0.890, 0.900, 0.910, 0.920,
         0.930, 0.940, 0.950, 0.960, 0.970, 0.980, 0.990, 1, 1.01, 1.02, 1.03, 1.04,
         1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17,
         1.18, 1.19, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.30,
         1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43,
         1.44, 1.449, 1.46, 1.471, 1.481, 1.493, 1.504, 1.515, 1.527, 1.538, 1.563,
         1.587, 1.613, 1.65, 1.68, 1.70, 1.73, 1.76, 1.80, 1.83, 1.84, 1.85, 1.855,
         1.86, 1.87, 1.89, 1.905, 1.923, 1.942, 1.961, 1.98, 2, 2.02, 2.041, 2.062,
         2.083, 2.105, 2.13, 2.15, 2.17, 2.19, 2.22, 2.24, 2.245, 2.25, 2.26, 2.27,
         2.29, 2.31, 2.33, 2.35, 2.37, 2.39, 2.41, 2.43, 2.46, 2.50, 2.52, 2.55, 2.565,
         2.58, 2.59, 2.60, 2.62, 2.675, 2.725, 2.778, 2.817, 2.833, 2.849, 2.865,
         2.882, 2.899, 2.915, 2.933, 2.95, 2.967, 2.985, 3.003, 3.021, 3.04, 3.058,
         3.077, 3.096, 3.115, 3.135, 3.155, 3.175, 3.195, 3.215, 3.236, 3.257, 3.279,
         3.30, 3.322, 3.345, 3.367, 3.39, 3.413, 3.436, 3.46, 3.484, 3.509, 3.534,
         3.559, 3.624, 3.732, 3.775, 3.847, 3.969, 4.099, 4.239, 4.348, 4.387, 4.444,
         4.505, 4.547, 4.56, 4.58, 4.719, 4.904, 5, 5.10]
    
    reals = [1.3801, 1.3509, 1.3339, 1.3249, 1.3203, 1.3194, 1.3185, 
         1.3177, 1.317, 1.3163, 1.3157, 1.3151, 1.3145, 1.314, 1.3135, 1.313, 
         1.3126, 1.3121, 1.3117, 1.3114, 1.311, 1.3106, 1.3103, 1.31, 1.3097, 
         1.3094, 1.3091, 1.3088, 1.3085, 1.3083, 1.308, 1.3078, 1.3076, 1.3073, 
         1.3071, 1.3069, 1.3067, 1.3065, 1.3062, 1.306, 1.3059, 1.3057, 1.3055, 
         1.3053, 1.3051, 1.3049, 1.3047, 1.3046, 1.3044, 1.3042, 1.304, 1.3039, 
         1.3037, 1.3035, 1.3033, 1.3032, 1.303, 1.3028, 1.3027, 1.3025, 1.3023, 
         1.3022, 1.302, 1.3019, 1.3017, 1.3015, 1.3014, 1.3012, 1.301, 1.3009, 
         1.3007, 1.3005, 1.3003, 1.3002, 1.3 , 1.2998, 1.2997, 1.2995, 1.2993, 
         1.2991, 1.299, 1.2988, 1.2986, 1.2984, 1.2982, 1.298, 1.2979, 1.2977, 
         1.2975, 1.2973, 1.2971, 1.2969, 1.2967, 1.2965, 1.2963, 1.2961, 1.2959, 
         1.2957, 1.2955, 1.2953, 1.2951, 1.2949, 1.2946, 1.2944, 1.2941, 1.2939, 
         1.2937, 1.2934, 1.2931, 1.2929, 1.2927, 1.2924, 1.2921, 1.292, 1.2918, 
         1.2916, 1.2914, 1.2912, 1.2909, 1.2903, 1.2897, 1.289, 1.2879, 1.287, 
         1.2863, 1.2853, 1.2843, 1.2828, 1.2816, 1.2811, 1.2807, 1.2805, 1.2802, 
         1.2797, 1.2788, 1.278, 1.2771, 1.2762, 1.2756, 1.275, 1.2744, 1.2736, 
         1.2728, 1.2718, 1.2707, 1.2694, 1.2677, 1.2663, 1.2648, 1.2633, 1.2609, 
         1.2591, 1.2587, 1.2582, 1.2573, 1.2564, 1.2545, 1.2525, 1.2504, 1.2482, 
         1.2459, 1.2435, 1.2409, 1.2382, 1.2337, 1.227, 1.2232, 1.2169, 1.2135, 
         1.2097, 1.2071, 1.2043, 1.1983, 1.1776, 1.1507, 1.1083, 1.0657, 1.0453, 
         1.0236, 1.0001, 0.9747, 0.9563, 0.9538, 0.9678, 0.9873, 1.0026, 1.018, 
         1.039, 1.0722, 1.1259, 1.2089, 1.3215, 1.4225, 1.4933, 1.5478, 1.597, 
         1.6336, 1.6477, 1.6405, 1.6248, 1.6108, 1.5905, 1.5714, 1.5559, 1.5396, 
         1.5241, 1.5086, 1.4949, 1.4827, 1.471, 1.4604, 1.4502, 1.4411, 1.4328, 
         1.4146, 1.3924, 1.385, 1.375, 1.3623, 1.3526, 1.3447, 1.3406, 1.3401, 
         1.3412, 1.3444, 1.3473, 1.3482, 1.3491, 1.347, 1.3379, 1.3325, 1.3268]
    
    imags = [2.0e-11,
         2.0e-11, 2.0e-11, 2.0e-11, 2.0e-11, 2.365e-11, 2.669e-11, 3.135e-11,
         4.14e-11, 6.268e-11, 9.239e-11, 1.325e-10, 1.956e-10, 2.861e-10,
         4.172e-10, 5.889e-10, 8.036e-10, 1.076e-09, 1.409e-09, 1.813e-09,
         2.289e-09, 2.839e-09, 3.461e-09, 4.159e-09, 4.93e-09, 5.73e-09,
         6.89e-09, 8.58e-09, 1.04e-08, 1.22e-08, 1.43e-08, 1.66e-08, 1.89e-08,
         2.09e-08, 2.4e-08, 2.9e-08, 3.44e-08, 4.03e-08, 4.3e-08, 4.92e-08,
         5.87e-08, 7.08e-08, 8.58e-08, 1.02e-07, 1.18e-07, 1.34e-07, 1.4e-07,
         1.43e-07, 1.45e-07, 1.51e-07, 1.83e-07, 2.15e-07, 2.65e-07, 3.35e-07,
         3.92e-07, 4.2e-07, 4.44e-07, 4.74e-07, 5.11e-07, 5.53e-07, 6.02e-07,
         7.55e-07, 9.26e-07, 1.12e-06, 1.33e-06, 1.62e-06, 2.0e-06, 2.25e-06,
         2.33e-06, 2.33e-06, 2.17e-06, 1.96e-06, 1.81e-06, 1.74e-06, 1.73e-06,
         1.7e-06, 1.76e-06, 1.82e-06, 2.04e-06, 2.25e-06, 2.29e-06, 3.04e-06,
         3.84e-06, 4.77e-06, 5.76e-06, 6.71e-06, 8.66e-06, 1.02e-05, 1.13e-05,
         1.22e-05, 1.29e-05, 1.32e-05, 1.35e-05, 1.33e-05, 1.32e-05, 1.32e-05,
         1.31e-05, 1.32e-05, 1.32e-05, 1.34e-05, 1.39e-05, 1.42e-05, 1.48e-05,
         1.58e-05, 1.74e-05, 1.98e-05, 3.442e-05, 5.959e-05, 0.0001028,
         0.0001516, 0.000203, 0.0002942, 0.0003987, 0.0004941, 0.0005532,
         0.0005373, 0.0005143, 0.0004908, 0.0004594, 0.0003858, 0.0003105,
         0.0002659, 0.0002361, 0.0002046, 0.0001875, 0.000165, 0.0001522,
         0.0001411, 0.0001302, 0.000131, 0.0001339, 0.0001377, 0.0001432,
         0.0001632, 0.0002566, 0.0004081, 0.000706, 0.001108, 0.001442,
         0.001614, 0.00164, 0.001566, 0.001458, 0.001267, 0.001023, 0.0007586,
         0.0005255, 0.0004025, 0.0003235, 0.0002707, 0.0002228, 0.0002037,
         0.0002026, 0.0002035, 0.0002078, 0.0002171, 0.0002538, 0.0003138,
         0.0003858, 0.0004591, 0.0005187, 0.0005605, 0.0005956, 0.0006259,
         0.000682, 0.000753, 0.0007685, 0.0007647, 0.0007473, 0.0007392,
         0.0007437, 0.0007543, 0.0008059, 0.001367, 0.003508, 0.01346,
         0.03245, 0.04572, 0.06287, 0.08548, 0.1198, 0.169, 0.221, 0.276, 0.312,
         0.347, 0.388, 0.438, 0.493, 0.554, 0.612, 0.625, 0.593, 0.539, 0.491, 0.438,
         0.372, 0.30, 0.238, 0.193, 0.158, 0.121, 0.103, 0.0836, 0.0668, 0.05312,
         0.04286, 0.03523, 0.02887, 0.02347, 0.01921, 0.01586, 0.01326, 0.0113,
         0.008146, 0.006672, 0.006966, 0.008248, 0.01112, 0.01471, 0.01867,
         0.02411, 0.02656, 0.0299, 0.03179, 0.0309, 0.03007, 0.02883, 0.0194,
         0.01347, 0.0124, 0.0122]
    
    f1 = interpolate.interp1d(wavelengths, reals)
    xnew = np.arange(0.305,5,0.01)
    reals_new = f1(xnew)   # use interpolation function returned by `interp1d`
    
    f2 = interpolate.interp1d(wavelengths, imags)
    imags_new = f2(xnew)
    
    reals = reals_new
    imags = imags_new
    wavelengths = xnew
    
    return reals, imags, wavelengths



def calc_optical_params(side_length,depth,reals,imags,wavelengths,plots=False,report_dims = False, ThreeBand = False):
    
    SSA_list = []
    Assy_list = []
    absXS_list = []
    MAC_list = []
    X_list = []
    
    V = 1.5*np.sqrt(3)*side_length**2*depth    # volume
    Area_total = 3 * side_length * (np.sqrt(3)*side_length+depth*2) #total surface area 
    Area = Area_total/4   # projected area
    apothem = (2*Area) / (depth*6) # apothem is distance from centre point to midpoint of a side for hexagon
    diameter = 2*apothem # midpoint of one side to midpoint of opposite side
                
    ar = depth/side_length
    delta = 0.3
    
    for i in np.arange(0,len(wavelengths),1):
      
    
        mr = reals[i]
        mi = imags[i]
        wl = wavelengths[i]
        
        #------------------------------------------------
        #---------- input tables (see Figs. 4 and 7) ----
        #------------------------------------------------
        # SSA parameterization
        a = [  0.457593 ,  20.9738 ] #for ar=1
        
        # SSA correction for AR != 1 (Table 2)
        nc1 = 3
        nc2 = 4
        c_ij = np.zeros(nc1*nc2*2).reshape((nc1,nc2,2))
        #   ---------- Plates ----------  
        c_ij[:,0,0] = [  0.000527060 ,  0.309748   , -2.58028  ]
        c_ij[:,1,0] = [  0.00867596  , -0.650188   , -1.34949  ]
        c_ij[:,2,0] = [  0.0382627   , -0.198214   , -0.674495 ]
        c_ij[:,3,0] = [  0.0108558   , -0.0356019  , -0.141318 ]
        #   --------- Columns ----------
        c_ij[:,0,1] = [  0.000125752 ,  0.387729   , -2.38400  ]
        c_ij[:,1,1] = [  0.00797282  ,  0.456133   ,  1.29446  ]
        c_ij[:,2,1] = [  0.00122800  , -0.137621   , -1.05868  ]
        c_ij[:,3,1] = [  0.000212673 ,  0.0364655  ,  0.339646 ]
        
        # diffraction g parameterization
        b_gdiffr = [ -0.822315 , -1.20125    ,  0.996653 ]
        
        # raytracing g parameterization ar=1
        p_a_eq_1 = [  0.780550 ,  0.00510997 , -0.0878268 ,  0.111549 , -0.282453 ]
        
        #---- g correction for AR != 1 (Also applied to AR=1 as plate) (Table 3)
        nq1 = 3
        nq2 = 7 
        q_ij = np.zeros(nq1*nq2*2).reshape((nq1,nq2,2))
        #   ---------- Plates ----------  
        q_ij[:,0,0] = [ -0.00133106  , -0.000782076 ,  0.00205422 ]
        q_ij[:,1,0] = [  0.0408343   , -0.00162734  ,  0.0240927  ]
        q_ij[:,2,0] = [  0.525289    ,  0.418336    , -0.818352   ]
        q_ij[:,3,0] = [  0.443151    ,  1.53726     , -2.40399    ]
        q_ij[:,4,0] = [  0.00852515  ,  1.88625     , -2.64651    ]
        q_ij[:,5,0] = [ -0.123100    ,  0.983854    , -1.29188    ]
        q_ij[:,6,0] = [ -0.0376917   ,  0.187708    , -0.235359   ]
        #   ---------- Columns ----------
        q_ij[:,0,1] = [ -0.00189096  ,  0.000637430 ,  0.00157383 ]
        q_ij[:,1,1] = [  0.00981029  ,  0.0409220   ,  0.00908004 ]
        q_ij[:,2,1] = [  0.732647    ,  0.0539796   , -0.665773   ]
        q_ij[:,3,1] = [ -1.59927     , -0.500870    ,  1.86375    ]
        q_ij[:,4,1] = [  1.54047     ,  0.692547    , -2.05390    ]
        q_ij[:,5,1] = [ -0.707187    , -0.374173    ,  1.01287    ]
        q_ij[:,6,1] = [  0.125276    ,  0.0721572   , -0.186466   ]
        
        #--------- refractive index correction of asymmetry parameter
        c_g = np.zeros(4).reshape(2,2)
        c_g[:,0] = [  0.96025050 ,  0.42918060 ]
        c_g[:,1] = [  0.94179149 , -0.21600979 ]
        #---- correction for absorption 
        s = [  1.00014  ,  0.666094 , -0.535922 , -11.7454 ,  72.3600 , -109.940 ]
        u = [ -0.213038 ,  0.204016 ]
        
        # -------- selector for plates or columns
        if ar > 1.:
            col_pla = 1 #columns
        else:
            col_pla = 0 #plates & compacts
            
        #------------------------------------------------
        #------------ Size parameters -------------------
        #------------------------------------------------
        
        #--- absorption size parameter (Fig. 4, box 1)
        Chi_abs = mi/wl*V/Area
        
        #----- scattering size parameter (Fig. 7, box 1)
        Chi_scat = 2.*np.pi*np.sqrt(Area/np.pi)/wl
        
        #------------------------------------------------
        #------------ SINGLE SCATTERING ALBEDO ----------
        #------------------------------------------------
        
        if Chi_abs > 0:
            w_1= 1.- a[0] * (1.-np.exp(-Chi_abs*a[1]))  #for AR=1 (Fig. 4, box 2)
            l=np.zeros(nc1)
            for i in range(nc2): l[:] += c_ij[:,i,col_pla] * np.log10(ar)**i  #(Fig. 4, box 3)
            D_w= l[0]*np.exp( -(np.log( Chi_abs )- l[2] )**2 / (2.*l[1]**2))/( Chi_abs *l[1]*np.sqrt(2.*np.pi)) #(Fig. 4, box 3)
            w = w_1 + D_w #(Fig. 4, box 4)
        else:
            w = 1.
        
        #------------------------------------------------
        #--------------- ASYMMETRY PARAMETER ------------
        #------------------------------------------------
        
        # diffraction g
        g_diffr = b_gdiffr[0] * np.exp(b_gdiffr[1]*np.log(Chi_scat)) + b_gdiffr[2] #(Fig. 7, box 2)
        g_diffr = max([g_diffr,0.5])
        
        # raytracing g at 862 nm
        g_1 = 0.
        for i in range(len(p_a_eq_1)): g_1 += p_a_eq_1[i]*delta**i #(Fig. 7, box 3)
        
        p_delta=np.zeros(nq1) 
        for i in range(nq2): p_delta += q_ij[:,i,col_pla]*np.log10(ar)**i #(Fig. 7, box 4)
        Dg = 0.
        for i in range(nq1): Dg += p_delta[i]*delta**i #(Fig. 7, box 4)
        g_rt = 2.*(g_1 + Dg)-1.  #(Fig. 7, box 5)
        
        #--------- refractive index correction of asymmetry parameter (Fig. 7, box 6)
        epsilon = c_g[0,col_pla]+c_g[1,col_pla]*np.log10(ar)
        mr1 = 1.3038 #reference value @ 862 nm band
        C_m = abs((mr1-epsilon)/(mr1+epsilon)*(mr+epsilon)/(mr-epsilon)) #abs function added according to corrigendum to the original paper
        
        #---- correction for absorption (Fig. 7, box 7)
        if Chi_abs > 0:
            C_w0 = 0.
            for i in range(len(s)): C_w0 += s[i]*(1.-w)**i
            k = np.log10(ar)*u[col_pla]
            C_w1 = k*w-k+1.    
            C_w = C_w0*C_w1
        else: C_w = 1.
        
        # raytracing g at required wavelength
        g_rt_corr = g_rt*C_m*C_w #(Fig. 7, box 9)
        
        #------ Calculate total asymmetry parameter and check g_tot <= 1 (Fig. 7, box 9)
        g_tot = 1./(2.*w)*( (2.*w-1.)*g_rt_corr + g_diffr )
        g_tot = min([g_tot,1.])
        
    
        absXS = Area*(1-((np.exp(-4*np.pi*mi*V))/(Area*wl)))
        MAC = absXS/V*914 # divide by volume*mass to give mass absorption coefficient
        
        SSA_list.append(w)
        Assy_list.append(g_tot)
        absXS_list.append(absXS)
        MAC_list.append(MAC)

    if ThreeBand:
        Band1MAC = np.mean(MAC_list[0:49])
        Band2MAC = np.mean(MAC_list[50:119])
        Band3MAC = np.mean(MAC_list[120:249])
        Band1SSA = np.mean(SSA_list[0:49])
        Band2SSA = np.mean(SSA_list[50:119])
        Band3SSA = np.mean(SSA_list[120:249])
        Band1Assy = np.mean(Assy_list[0:49])
        Band2Assy = np.mean(Assy_list[50:119])
        Band3Assy = np.mean(Assy_list[120:249])
        Band1absXS = np.mean(absXS_list[0:49])
        Band2absXS = np.mean(absXS_list[50:119])
        Band3absXS = np.mean(absXS_list[120:249])

        ThreeBandAssy = [Band1Assy, Band2Assy, Band3Assy]
        ThreeBandSSA = [Band1SSA, Band2SSA, Band3SSA]
        ThreeBandabsXS = [Band1absXS, Band2absXS, Band3absXS]
        ThreeBandMAC = [Band1MAC,Band2MAC,Band3MAC]

    else:
        ThreeBandSSA = 0
        ThreeBandMAC = 0
        ThreeBandAssy = 0
        ThreeBandabsXS = 0

    if plots:
        plt.figure(1)    
        plt.plot(wavelengths,SSA_list),plt.ylabel('SSA'),plt.xlabel('Wavelength (um)'),plt.grid(b=None)
        plt.figure(2)
        plt.plot(wavelengths,Assy_list),plt.ylabel('Assymetry Parameter'),plt.xlabel('Wavelength (um)'),plt.grid(b=None)
        plt.figure(3)
        plt.plot(wavelengths,MAC_list),plt.ylabel('Mass Absorption Cross Section'),plt.xlabel('Wavelength (um)'),plt.grid(b=None)
    
    if report_dims:
        print('Width of hexagonal plane = ',np.round(diameter/10000,2),' (cm)')
        print('depth of hexagonal column = ',depth/10000,' (cm)')
        print('aspect ratio = ',ar)
        print('ice crystal volume = ', np.round(V*1e-12,2),' (cm^3)')
    
    return Assy_list,SSA_list,MAC_list,depth,side_length, ThreeBandSSA, ThreeBandAssy, ThreeBandabsXS, \
           ThreeBandMAC


def net_cdf_updater(filepath,filepath3Band, Assy_list, SSA_list, MAC_list, depth, side_length, density, ThreeBandSSA,
                    ThreeBandAssy, ThreeBandMAC, ThreeBand = False):

    if ThreeBand:
        filepathIN = filepath3Band
        MAC_IN = ThreeBandMAC
        SSA_IN = ThreeBandSSA
        Assy_IN = ThreeBandAssy

    else:
        filepathIN = filepath
        MAC_IN = np.squeeze(MAC_list)
        SSA_IN = np.squeeze(SSA_list)
        Assy_IN = np.squeeze(Assy_list)


    icefile = pd.DataFrame()
    icefile['asm_prm'] = Assy_IN
    icefile['ss_alb'] = SSA_IN
    icefile['ext_cff_mss'] = MAC_IN
    icefile = icefile.to_xarray()
    icefile.attrs['medium_type'] = 'air'
    icefile.attrs['description'] = 'Optical properties for ice grain: hexagonal column of side length {}um and length {}um'.format(
        str(side_length), str(depth))
    icefile.attrs['psd'] = 'monodisperse'
    icefile.attrs['side_length_um'] = depth
    icefile.attrs['density_kg_m3'] = density
    icefile.attrs[
        'origin'] = 'Optical properties derived from geometrical optics calculations'
    icefile.to_netcdf(str(filepathIN + 'ice_wrn_{}_{}_3Band.nc'.format(str(side_length), str(depth))))

    return 

###############################################################################
##########################  FUNCTON CALLS ####################################

reals,imags,wavelengths = preprocess_RI()

for side_length in np.arange(1000,15000,500):
    for depth in np.arange(1000,15000,500):
        Assy_list, SSA_list, MAC_list, depth, side_length, ThreeBandSSA, ThreeBandAssy, ThreeBandabsXS, \
        ThreeBandMAC = calc_optical_params(side_length,depth, reals,imags,wavelengths,plots=False,
                                                           report_dims = False, ThreeBand = True)
        net_cdf_updater(filepath,filepath3Band,Assy_list,SSA_list,MAC_list, depth, side_length, 917,
                       ThreeBandSSA, ThreeBandAssy, ThreeBandMAC, ThreeBand=True)
