
# def snicar8d_mie(BND_TYP, DIRECT, CLEAR, CLOUDY, APRX_TYP,DELTA,coszen,R_sfc,dz,rho_snw,nbr_aer,MSSsoot1,MSSsoot2,MSSdust1,MSSdust2,
#                  MSSdust3,MSSdust4,MSS_ash,MSSGRISdust,MSSsnowalg,MSSglacieralg1,MSSglacieralg2,MSSglacieralg3,
#                  MSSglacieralg4,MSSglacieralg5,FILEsoot1, FILESsoot2, FILEash, FILEGRISdust FILEdust1,FILEdust2,FILEdust3,FILEdust4,FILEsnowalg,FILEglacieralg1,
#                  FILEglacieralg2,FILEglacieralg3,FILEglacieralg4,FILEglacieralg5)


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

""" temporary variable assignment for testing: will be defined in driver function eventually """

BND_TYP = 1
coszen = 0.57
DIRECT = True
DELTA = False
APRX_TYP = 1
R_sfc = 0.15
rds_snw = np.array([500,500,500,500,500])
rho_snw = np.array([300,300,300,300,300])
nbr_lyr = 5
nbr_aer = 14
dz = np.array([0.001, 0.01, 0.01, 0.01, 0.01])

FILEsoot1 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/miecot_slfsot_ChC90_dns_1317.nc")
FILEsoot2 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/mie_sot_ChC90_dns_1317.nc")
FILEash = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/volc_ash_mtsthelens_20081011.nc")
FILEGRISdust = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/GRISdust_PSD.nc")
FILEdust1 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/aer_dst_bln_20060904_01.nc")
FILEdust2 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/aer_dst_bln_20060904_02.nc")
FILEdust3 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/aer_dst_bln_20060904_03.nc")
FILEdust4 = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Mie_files/aer_dst_bln_20060904_04.nc")
FILEsnowalg = xr.open_dataset("/home/joe/Code/BioSNICAR_GO/Algal_Optical_Props/snw_alg_1.nc")
FILEglacieralg1 = xr.open_dataset("/home/joe/Desktop/BioSNICAR_GO_CW/Algal_Optical_Props/CW_5_algae_geom_5_40.nc")
FILEglacieralg2 = xr.open_dataset("/home/joe/Desktop/BioSNICAR_GO_CW/Algal_Optical_Props/CW_5_algae_geom_5_40.nc")
FILEglacieralg3 = xr.open_dataset("/home/joe/Desktop/BioSNICAR_GO_CW/Algal_Optical_Props/CW_5_algae_geom_5_40.nc")
FILEglacieralg4 = xr.open_dataset("/home/joe/Desktop/BioSNICAR_GO_CW/Algal_Optical_Props/CW_5_algae_geom_5_40.nc")
FILEglacieralg5 = xr.open_dataset("/home/joe/Desktop/BioSNICAR_GO_CW/Algal_Optical_Props/CW_5_algae_geom_5_40.nc")

MSSsoot1 = np.array([0,0,0,0,0])
MSSsoot2 = np.array([0,0,0,0,0])
MSSash = np.array([0,0,0,0,0])
MSSGRISdust = np.array([0,0,0,0,0])
MSSdust1 = np.array([0,0,0,0,0])
MSSdust2 = np.array([0,0,0,0,0])
MSSdust2 = np.array([0,0,0,0,0])
MSSdust3 = np.array([0,0,0,0,0])
MSSdust4 = np.array([0,0,0,0,0])
MSSsnowalg = np.array([0,0,0,0,0])
MSSglacieralg1 = np.array([0,0,0,0,0])
MSSglacieralg2 = np.array([0,0,0,0,0])
MSSglacieralg3 = np.array([0,0,0,0,0])
MSSglacieralg4 = np.array([0,0,0,0,0])
MSSglacieralg5 = np.array([0,0,0,0,0])



""" BEGIN FUNCTION """

# set working directory (location of netcdf library)
dir_base = "/home/joe/Code/BioSNICAR_GO/"
dir_alg = "Algal_Optical_Props/"
dir_mie_files = "Mie_files/"

# retrieve wavelength from arbitrary choice of netcdf file
temp = xr.open_dataset(str(dir_base+dir_mie_files+"ice_wrn_0200.nc"))
wvl = np.array(temp['wvl'].values)
wvl = wvl*1e6
nbr_wvl = len(wvl)

# set reflectance of underlying surface
R_sfc = [R_sfc for _ in range(nbr_wvl)]
R_sfc=np.array(R_sfc)

# Incoming Irradiance
# calculate mu_not
#mu_not = np.cos((slr_znt / 360) * 2 * np.pi() # convert radians if required
mu_not = coszen
flx_slr = []

if DIRECT:

    with open(str(dir_base + dir_mie_files + "mlw_sfc_flx_frc_clr.txt")) as file:
        for line in file:
            line = float(line.rstrip("\n"))
            flx_slr.append(line)
    flx_slr = np.array(flx_slr)
    flx_slr[flx_slr==0]=1e-30
    Fs = flx_slr / (mu_not * np.pi)

    Fd = np.zeros(nbr_wvl)

else:
    with open(str(dir_base + dir_mie_files + "mlw_sfc_flx_frc_cld.txt")) as file:
        for line in file:
            line = float(line.rstrip("\n"))
            flx_slr.append(line)
    flx_slr = np.array(flx_slr)
    flx_slr[flx_slr==0]=1e-30

    Fd = [flx_slr[i]/mu_not*np.pi for i in range(nbr_wvl)]
    Fs = np.zeros(nbr_wvl)

# Read in ice optical properties
# set string stubs for reading in ice optical files

fl_stb1 = "ice_wrn_"
fl_stb2 = ".nc"

#set up empty arrays
SSA_snw = np.empty([nbr_lyr, nbr_wvl])
MAC_snw = np.empty([nbr_lyr, nbr_wvl])
g_snw = np.empty([nbr_lyr, nbr_wvl])

for i in np.arange(0,nbr_lyr,1):

    if rds_snw[i] ==0:
        print("ERROR: ICE GRAIN RADIUS SET TO ZERO")

    else:
        s1 = "{}".format(str(rds_snw[i]))
        s1 = s1.rjust(4,'0')
        FILE_ice = str(dir_base+dir_mie_files+fl_stb1+s1+fl_stb2)

# read in single scattering albedo, MAC and g for ice crystals in each layer
    temp = xr.open_dataset(FILE_ice)
    SSA = temp['ss_alb'].values
    SSA_snw[i,:] = SSA

    ext_cff_mss = temp['ext_cff_mss'].values
    MAC_snw[i,:] = ext_cff_mss

    asm_prm = temp['asm_prm'].values
    g_snw[i,:] = asm_prm

# read in aerosol optical properties

SSAaer = np.zeros([nbr_aer,nbr_wvl])

SSAaer[0,:] = FILEsoot1['ss_alb'].values
SSAaer[1,:] = FILEsoot2['ss_alb'].values
SSAaer[2,:] = FILEash['ss_alb'].values
SSAaer[3,:] = FILEdust1['ss_alb'].values
SSAaer[4,:] = FILEdust2['ss_alb'].values
SSAaer[5,:] = FILEdust3['ss_alb'].values
SSAaer[6,:] = FILEdust4['ss_alb'].values
SSAaer[7,:] = FILEGRISdust['ss_alb'].values
SSAaer[8,:] = FILEsnowalg ['ss_alb'].values
SSAaer[9,:] = FILEglacieralg1['ss_alb'].values
SSAaer[10,:] = FILEglacieralg2['ss_alb'].values
SSAaer[11,:] = FILEglacieralg3['ss_alb'].values
SSAaer[12,:] = FILEglacieralg4['ss_alb'].values
SSAaer[13,:] = FILEglacieralg5['ss_alb'].values

MACaer = np.zeros([nbr_aer, nbr_wvl])

MACaer[0,:] = FILEsoot1['ext_cff_mss'].values
MACaer[1,:] = FILEsoot2['ext_cff_mss'].values
MACaer[2,:] = FILEash['ext_cff_mss'].values
MACaer[3,:] = FILEdust1['ext_cff_mss'].values
MACaer[4,:] = FILEdust2['ext_cff_mss'].values
MACaer[5,:] = FILEdust3['ext_cff_mss'].values
MACaer[6,:] = FILEdust4['ext_cff_mss'].values
MACaer[7,:] = FILEGRISdust['ext_cff_mss'].values
MACaer[8,:] = FILEsnowalg['ext_cff_mss'].values
MACaer[9,:] = FILEglacieralg1['ext_cff_mss'].values
MACaer[10,:] = FILEglacieralg2['ext_cff_mss'].values
MACaer[11,:] = FILEglacieralg3['ext_cff_mss'].values
MACaer[12,:] = FILEglacieralg4['ext_cff_mss'].values
MACaer[13,:] = FILEglacieralg5['ext_cff_mss'].values

Gaer = np.zeros([nbr_aer,nbr_wvl])

Gaer[0,:] = FILEsoot1['asm_prm'].values
Gaer[1,:] = FILEsoot2['asm_prm'].values
Gaer[2,:] = FILEash['asm_prm'].values
Gaer[3,:] = FILEdust1['asm_prm'].values
Gaer[4,:] = FILEdust2['asm_prm'].values
Gaer[5,:] = FILEdust3['asm_prm'].values
Gaer[6,:] = FILEdust4['asm_prm'].values
Gaer[7,:] = FILEGRISdust['asm_prm'].values
Gaer[8,:] = FILEsnowalg['asm_prm'].values
Gaer[9,:] = FILEglacieralg1['asm_prm'].values
Gaer[10,:] = FILEglacieralg2['asm_prm'].values
Gaer[11:] = FILEglacieralg3['asm_prm'].values
Gaer[12,:] = FILEglacieralg4['asm_prm'].values
Gaer[12,:] = FILEglacieralg5['asm_prm'].values

# load mass concentrations per layer into numpy array (one row per layer, one column per umpurity)
# and convert to kg/kg unit

MSSaer = np.zeros([nbr_lyr, nbr_aer])
MSSaer[0:nbr_lyr,0] = MSSsoot1
MSSaer[0:nbr_lyr,1] = MSSsoot2
MSSaer[0:nbr_lyr,2] = MSSash
MSSaer[0:nbr_lyr,3] = MSSdust1
MSSaer[0:nbr_lyr,4] = MSSdust2
MSSaer[0:nbr_lyr,5] = MSSdust3
MSSaer[0:nbr_lyr,6] = MSSdust4
MSSaer[0:nbr_lyr,7] = MSSGRISdust
MSSaer[0:nbr_lyr,8] = MSSsnowalg
MSSaer[0:nbr_lyr,9] = MSSglacieralg1
MSSaer[0:nbr_lyr,10] = MSSglacieralg2
MSSaer[0:nbr_lyr,11] = MSSglacieralg3
MSSaer[0:nbr_lyr,12] = MSSglacieralg4
MSSaer[0:nbr_lyr,13] = MSSglacieralg5
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
direct = np.zeros([nbr_lyr,nbr_wvl])
F_net = np.zeros([nbr_lyr, nbr_wvl])
F_btm_net = np.zeros([1,nbr_wvl])
F_top_net = np.zeros([1,nbr_wvl])
intensity = np.zeros([nbr_lyr, nbr_wvl])
F_top_pls = np.zeros([1,nbr_wvl])
F_up = np.zeros([nbr_lyr,nbr_wvl])
F_down = np.zeros([nbr_lyr,nbr_wvl])
F_net2 = np.zeros([nbr_lyr, nbr_wvl])
intensity2 = np.zeros([nbr_lyr, nbr_wvl])
intensity2_top = np.zeros(nbr_wvl)
F_abs = np.zeros([nbr_lyr, nbr_wvl])


# for each layer, the layer mass (L) is density * layer thickness
# for each layer the optical depth is the layer mass * the mass extinction coefficient
# first for the ice in each layer
for i in range(nbr_lyr):
    L_snw[i] = rho_snw[i] * dz[i]
    tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]

# then for the LAPs in each layer
for i in range(nbr_lyr):
    for j in range(nbr_aer):
        L_aer[i, j, :] = np.multiply(L_snw[i], MSSaer[i, j])
        tau_aer[i, j, :] = np.multiply(L_aer[i, j, :], MACaer[j, :])

        tau_sum = tau_sum + tau_aer[i, j, :]
        SSA_sum = SSA_sum + (tau_aer[i, j, :] * SSAaer[j, :])
        g_sum = g_sum + (tau_aer[i, j, :] * SSAaer[j, :] * Gaer[j, :])

# finally, for each layer calculate the effective SSA, tau and g for the snow+LAP
for i in range(nbr_lyr):
    tau[i,:] = tau_sum[i,:] + tau_snw[i,:]
    SSA[i,:] = (1/tau[i,:]) * (SSA_sum[i,:] + SSA_snw[i,:] * tau_snw[i,:])
    g[i, :] = (1 / (tau[i, :] * (SSA[i, :]))) * (g_sum[i,:] + (g_snw[i, :] * SSA_snw[i, :] * tau_snw[i, :]))


############################################
# PERFORM DELTA TRANSFORMATION IF REQUIRED
############################################
# The star represents the delta transformed quantity
# if no delat transformation is applied, the starred quantity
# is equal to the unstarred quantity

if DELTA:
    g_star = g/(1+g)
    SSA_star = ((1-(g**2))*SSA)/(1-(SSA*(g**2)))
    tau_star = (1-(SSA*(g**2)))*tau

else:
    g_star = g
    SSA_star = SSA
    tau_star = tau


# CALCULATE TOTAL OPTICAL DEPTH OF ENTIRE COLUMN
# i.e. tau_clm = total optical depth from upper boundary 
# to upper boundary of layer n. This is therefore a cumulative
# quantity - subsequently lower layers contain the sum of the
# # optical depth of all overlying layers

tau_clm = np.zeros([nbr_lyr,nbr_wvl])
for i in np.arange(1,nbr_lyr,1):
    #start loop from 2nd layer, i.e. index = 1
    tau_clm[i,:] = tau_clm[i-1,:]+tau_star[i-1,:]

# SET BOUNDARY CONDITION: BOTTOM BOUNDARY
# calculate radiation reflected skywards by underlying surface (i.e. lower model boundary)
# remainder is lost

S_sfc = R_sfc * mu_not * np.exp(-(tau_clm[nbr_lyr-1,:] + tau_star[nbr_lyr-1,:])/mu_not)*np.pi * Fs


######################################################
# Apply Two-Stream Approximation (Toon et al, table 1)
######################################################
"""
Three 2-stream approximations are available: Eddington,
Quadrature and hemispheric mean. The equations for each
approximation are provided in Toon et al. (1989) Table 1.

The hemispheric mean scheme is derived by assuming that the
phase function is equal to 1  + g  in the forward scattering 
hemisphere and to 1  - g  in the backward scattering hemisphere. 
The asymmetry parameter is g. The hemispheric mean is only
useful for infrared wavelengths

"""

if APRX_TYP == 1:
    #apply Eddington approximation
    gamma1 = (7-(SSA_star * (4+(3*g_star))))/4
    gamma2 = -(1-(SSA_star*(4-(3*g_star))))/4
    gamma3 = (2-(3*g_star*mu_not))/4
    gamma4 = 1-gamma3
    mu_one = 0.5

elif APRX_TYP==2:
    #apply quadrature approximation
    gamma1 = np.sqrt(3)*(2-(SSA_star*(1+g_star)))/2
    gamma2 = SSA_star * np.sqrt(3)*(1-g_star)/2
    gamma3 = (1-(np.sqrt(3)*g_star*mu_not))/2
    gamma4 = 1-gamma3
    mu_one = 1/np.sqrt(3)

elif APRX_TYP==3:
    #apply hemispheric mean approximation
    gamma1 = 2 - (SSA_star*(1+g_star))
    gamma2 = SSA_star*(1-g_star)
    gamma3 = (1-(np.sqrt(3) * g_star*mu_not))/2
    gamma4 = 1-gamma3
    mu_one = 0.5


# Toon et al equation 21 and 22
# Note that the values of lam and GAMMA depend upon gamma1 and gamma2, which
# vary depending upon the two-stream approximation used
# variable "lambda" renamed "lam" to avoid confusion with lambda function
lam = np.sqrt(abs((gamma1**2)-(gamma2**2)))
GAMMA = gamma2/(gamma1+lam)

# calculate coefficients required for tridiagonal matrix calculation
# (Toon et al Equation 44)
e1 = 1+(GAMMA*np.exp(-lam*tau_star))
e2 = 1-(GAMMA*np.exp(-lam*tau_star))
e3 = GAMMA+np.exp(-lam*tau_star)
e4 = GAMMA-np.exp(-lam*tau_star)


######################################
# Calculate C-functions
######################################

# C is calculated at the top and bottom of each layer, i,
# see Toon equations 23 and 24

""" N.B. consider adding in stability check here as per Flanner's Matlab code """

C_pls_btm = np.zeros([nbr_lyr,nbr_wvl])
C_mns_btm = np.zeros([nbr_lyr,nbr_wvl])
C_pls_top = np.zeros([nbr_lyr, nbr_wvl])
C_mns_top = np.zeros([nbr_lyr, nbr_wvl])

for i in range(nbr_lyr):

    if np.sum(Fs) > 0.0:

        C_pls_btm[i,:] = (SSA_star[i,:]*np.pi*Fs*np.exp(-(tau_clm[i,:]+tau_star[i,:])/mu_not)*
                            (((gamma1[i,:]-(1/mu_not))*gamma3[i,:])+(gamma4[i,:]*gamma2[i,:])))\
                            /((lam[i,:]**2)-(1/(mu_not**2)))

        C_mns_btm[i,:] = (SSA_star[i,:]*np.pi*Fs*
                            np.exp(-(tau_clm[i,:]+tau_star[i,:])/mu_not) * (((gamma1[i,:]+(1/mu_not))*gamma4[i,:])+
                            (gamma2[i,:]*gamma3[i,:])))/((lam[i,:]**2)-(1/mu_not**2))

        C_pls_top[i,:] = (SSA_star[i,:] * np.pi * Fs * np.exp(-tau_clm[i,:]/mu_not)* ((gamma1[i,:] - (1/mu_not))
                            * gamma3[i,:] + (gamma4[i,:]*gamma2[i,:])))/((lam[1,:]**2)-(1/mu_not**2))

        C_mns_top[i,:] = (SSA_star[i,:] * np.pi * Fs * np.exp(-tau_clm[i,:]/mu_not) * ((gamma1[i,:]+(1/mu_not))
                            * gamma4[i,:] + (gamma2[i,:]*gamma3[i,:])))/((lam[i,:]**2)-(1/mu_not**2))

    else:
        # no direct-beam flux:
        C_pls_btm[i,:] = 0
        C_mns_btm[i,:] = 0
        C_pls_top[i,:] = 0
        C_mns_top[i,:] = 0


# Toon equations 41-43.
# Boundary values for i=1 and i=2nbr_lyr, specifics for i=odd and i=even
# Set up lists
A = np.zeros([2*nbr_lyr,nbr_wvl])
B = np.zeros([2*nbr_lyr,nbr_wvl])
D = np.zeros([2*nbr_lyr,nbr_wvl])
E = np.zeros([2*nbr_lyr,nbr_wvl])


###########################################
# Initialize tridiagonal matrix solution
###########################################

# expanding the number of layers to 2*nbr_lyr

####### NB HARDCODING LAYER NUMBERS IN THE FOLLOWING A,B,D,E CALCULATIONS MAKES THE CODE 
####### APPLICABLE TO 5 LAYER CONFIGURATIONS ONLY!!
# TODO: Fix this!

for i in np.arange(0,2*nbr_lyr,1):

    #TOP LAYER    
    if i==0:
        A[0,:] = 0.0
        B[0,:] = e1[0,:]
        D[0,:] = -e2[0,:]
        E[0,:] = Fd-C_mns_top[0,:]

    # BOTTOM LAYER
    elif i== 2*nbr_lyr-1:
        A[i,:] = e1[nbr_lyr-1,:]-(R_sfc * e3[nbr_lyr-1,:])
        B[i,:] = e2[nbr_lyr-1,:]-(R_sfc * e4[nbr_lyr-1,:])
        D[i,:] = 0.0
        E[i,:] = S_sfc[:] - C_pls_btm[nbr_lyr-1,:] + (R_sfc * C_mns_btm[nbr_lyr-1,:])
  
    # LAYERS 2,4
    elif (i == 2)|(1==4):
        n = int(i/2)
        A[i,:] = (e2[n,:] * e3[n,:])-(e4[n,:] * e1[n,:])
        B[i,:] = (e1[n,:] * e1[n+1,:])-(e3[n,:] * e3[n+1,:])
        D[i,:] = (e3[n,:] * e4[n+1,:])-(e1[n,:] * e2[n+1,:])
        E[i,:] = (e3[n,:] * (C_pls_top[n+1,:] - C_pls_btm[n,:])) +  (e1[n,:] * (C_mns_btm[n,:] - C_mns_top[n+1,:]))

    # LAYERS 1,3
    elif (i == 1)|(1==3):
        n = int(np.floor(i/2))
        A[i,:] = (e2[n+1,:] * e1[n,:])-(e3[n,:] * e4[n+1,:])
        B[i,:] = (e2[n,:] * e2[n+1,:])-(e4[n,:] * e4[n+1,:])
        D[i,:] = (e1[n+1,:] * e4[n+1,:])-(e2[n+1,:] * e3[n+1,:])
        E[i,:] = (e2[n+1,:] * (C_pls_top[n+1,:] - C_pls_btm[n,:])) + (e4[n+1,:] * (C_mns_top[n+1,:] - C_mns_btm[n,:]))


# Now the actual tridiagonal matrix solving. Simply dividing A/B and E/B 
# throws an exception due to division by zero. Here we use numpy's nan_to_num
# function to achieve the division where possible and replace nans with zeros.
# We also set numpy to ignore the division error.

# for bottom layer only
# Toon et al Eq 45
AS = np.zeros([2*nbr_lyr,nbr_wvl])
DS = np.zeros([2*nbr_lyr,nbr_wvl])

np.seterr(divide='ignore',invalid='ignore')
AS[2*nbr_lyr-1,:] = np.nan_to_num(A[2*nbr_lyr-1,:]/B[2*nbr_lyr-1,:])
DS[2*nbr_lyr-1,:] = np.nan_to_num(E[2*nbr_lyr-1,:]/B[2*nbr_lyr-1,:])

# for all layers above bottom layer, starting at second-to-bottom and progressing towards
# surface:
# Toon et al Eq 46
X = np.zeros([nbr_lyr*2-1,nbr_wvl])
for i in np.arange(2*nbr_lyr-2,-1, -1):
    X[i,:] = 1/(B[i,:]-(D[i,:] * AS[i+1,:]))
    AS[i,:] = np.nan_to_num(A[i,:]*X[i,:])
    DS[i,:] = np.nan_to_num((E[i,:]-(D[i,:]*DS[i+1,:]))*X[i,:])

# then for all layers, progressing from surface to bottom
# Toon et al Eq 47
Y = np.zeros([nbr_lyr*2,nbr_wvl])

for i in np.arange(0,2*nbr_lyr,1):
    if i ==0:
        Y[0,:] = DS[0,:]
    else:
        Y[i,:] = DS[i,:] - (AS[i,:]*Y[i-1,:])
   

#############################################################
# CALCULATE DIRECT BEAM FLUX AT BOTTOM OF EACH LAYER


# loop through layers
for i in np.arange(0,nbr_lyr,1):
    # (Toon et al. eq 50)
    direct[i,:] = mu_not * np.pi * Fs * np.exp(-(tau_clm[i,:] + tau_star[i,:]) / mu_not)

    # net flux (positive upward = F_up - F_down) at the base of each layer (Toon et al. Eq 48)
    F_net[i,:] = (Y[2*i-1,:] * (e1[i,:] - e3[i,:])) + (Y[2*i-1] * (e2[i,:] - e4[i,:])) + C_pls_btm[i,:] - C_mns_btm[i,:] - direct[i,:]

    # mean intensity at the base of each layer (Toon et al. Eq 49)
    intensity[i,:] = (1/mu_one) * (Y[2*i,:] * (e1[i,:] + e3[i,:]) + Y[2*i-1,:] * (e2[i,:] + e4[i,:]) + C_pls_btm[i,:] + C_mns_btm[i,:]) + (direct[i,:]/mu_not)
    intensity[i, :] = intensity[i, :] / (4 * np.pi)


# Upward flux at upper model boundary (Toon et al Eq 31)
F_top_pls = (Y[0,:] * (np.exp(-lam[0,:] * tau_star[0,:]) + GAMMA[0,:])) + (Y[1,:] * (np.exp(-lam[0,:] * tau_star[0,:])-GAMMA[0,:])) + C_pls_top[0,:]


for i in np.arange(0,nbr_lyr,1):
    # Upward flux at the bottom of each layer interface (Toon et al. Eq31)
    F_up[i,:] = Y[2*i,:] * (np.exp(0) + GAMMA[i,:] * np.exp(-lam[i,:] * tau_star[i,:])) + Y[2*i-1,:] * (np.exp(0) - GAMMA[i,:] * np.exp(-lam[i,:] * tau_star[i,:])) + C_pls_btm[i,:]

    # Downward flux at the bottom of each layer interface (Toon et al. Eq32) plus direct beam component
    F_down[i,:] = Y[2*i,:] * (GAMMA[i,:] * np.exp(0) + np.exp(-lam[i,:] * tau_star[i,:])) + Y[2*i-1,:] * (GAMMA[i,:] * np.exp(0) - np.exp(-lam[i,:] * tau_star[i,:])) + C_mns_btm[i,:] + direct[i,:];

    # Derived net (upward-downward) flux (should equal F_net)
    F_net2[i,:] = F_up[i,:] - F_down[i,:]

    # surface planar intensity
    intensity2_top[:] = F_top_pls + ((mu_not * np.pi * Fs) + Fd)

# Net flux at lower model boundary = bulk transmission through entire media
# = energy absorbed by underlying surface
F_btm_net[0,:] = -F_net[-1,:]

# Hemispheric wavelength-dependent albedo
albedo = F_top_pls/ ((mu_not * np.pi * Fs)+Fd)

# Net flux at upper model boundary
F_top_net[0,:] = F_top_pls - ((mu_not * np.pi * Fs) + Fd)

# absorbed flux in each layer (negative if there is net emission (bnd_typ = 4))
for i in np.arange(0,nbr_lyr,1):
    if i ==1:
        F_abs[0,:] = F_net[0,:]-F_top_net
    else:
        F_abs[i,:] = F_net[i,:] - F_net[i-1,:]

# set indices for constraining calculations to VIS and NIR bands
if BND_TYP == 1:
    vis_max_idx = 39
    nir_max_idx = len(wvl)

elif (BND_TYP == 2) | (BND_TYP == 3):
    vis_max_idx = 0
    nir_max_idx = len(wvl)

# Spectrally-integrated absorption in each layer:
F_abs_slr = np.sum(F_abs,axis=1)

F_abs_vis = np.zeros([nbr_lyr, 1])
F_abs_nir = np.zeros([nbr_lyr, 1])

for i in np.arange(0,nbr_lyr,1):
    F_abs_vis[i] = np.sum(F_abs[i,0:vis_max_idx])
    F_abs_nir[i] = np.sum(F_abs[i,vis_max_idx:nir_max_idx])


# Spectrally - integrated absorption by underlying surface:
F_abs_btm = np.sum(F_btm_net,axis=1)
F_abs_vis_btm = np.sum(F_btm_net[0:vis_max_idx],axis=1)
F_abs_nir_btm = np.sum(F_btm_net[vis_max_idx:nir_max_idx],axis=1)

# Calculate radiative heating rate in kelvin per second.
# Multiply by 3600 to convert to K per hour
# specfic heta capacity of ice = 2117 J kg-1 K-1
heat_rt = F_abs_slr / (L_snw * 2117) # [K / s]
heat_rt = heat_rt * 3600 # [K / hr]


# Energy conservation check:
# % Incident direct + diffuse radiation equals(absorbed + transmitted + bulk_reflected)
energy_sum = (mu_not * np.pi * Fs) + Fd - (np.sum(F_abs,axis=0) + F_btm_net + F_top_pls)

# spectrally-integrated terms:
# energy conservation total error
energy_error = abs(np.sum(energy_sum))

if energy_error > 1e-10:
    energy_conservation_error = np.sum(abs(energy_sum))
    print(f"CONSERVATION OF ENERGY ERROR OF {energy_conservation_error}")

# total incident insolation(Wm - 2)
total_insolation = np.sum((mu_not * np.pi * Fs) + Fd)

# energy absorbed by all snow layers
energy_abs_all_lyrs = np.sum(np.sum(F_abs))

# energy absorbed by underlying substrate
energy_abs_under_sfc = np.sum(F_btm_net)

# Spectrally - integrated solar, visible, and NIR albedos:
alb = np.sum(flx_slr * albedo) / np.sum(flx_slr)


#alb_vis = sum(flx_slr(1:vis_max_idx).*albedo(1:vis_max_idx))/...
      #    sum(flx_slr(1:vis_max_idx));

# alb_nir = sum(flx_slr(vis_max_idx + 1:nir_max_idx).*albedo(vis_max_idx + 1: nir_max_idx)) / ...
# sum(flx_slr(vis_max_idx + 1: nir_max_idx));

# % Diagnostic
# for comparing 470 - band solutions with 5-band solutions:
#     %
# Spectrally - integrated
# 5 - band
# albedo:
# if (BND_TYP == 1)
# bnd1a=1;
# bnd1b=40;
# bnd2a=41;
# bnd2b=70;
# bnd3a=71;
# bnd3b=90;
# bnd4a=91;
# bnd4b=120;
# bnd5a=121;
# bnd5b=470;

# bnd6a=1;
# bnd6b=40;
# bnd7a=41;
# bnd7b=89;
# bnd8a=90;
# bnd8b=470;

# alb1 = sum(flx_slr(bnd1a:bnd1b).*albedo(bnd1a: bnd1b)) / sum(flx_slr(bnd1a: bnd1b));
# alb2 = sum(flx_slr(bnd2a:bnd2b).*albedo(bnd2a: bnd2b)) / sum(flx_slr(bnd2a: bnd2b));
# alb3 = sum(flx_slr(bnd3a:bnd3b).*albedo(bnd3a: bnd3b)) / sum(flx_slr(bnd3a: bnd3b));
# alb4 = sum(flx_slr(bnd4a:bnd4b).*albedo(bnd4a: bnd4b)) / sum(flx_slr(bnd4a: bnd4b));
# alb5 = sum(flx_slr(bnd5a:bnd5b).*albedo(bnd5a: bnd5b)) / sum(flx_slr(bnd5a: bnd5b));
# alb6 = sum(flx_slr(bnd6a:bnd6b).*albedo(bnd6a: bnd6b)) / sum(flx_slr(bnd6a: bnd6b));
# alb7 = sum(flx_slr(bnd7a:bnd7b).*albedo(bnd7a: bnd7b)) / sum(flx_slr(bnd7a: bnd7b));
# alb8 = sum(flx_slr(bnd8a:bnd8b).*albedo(bnd8a: bnd8b)) / sum(flx_slr(bnd8a: bnd8b));
# end;

# % Spectrally - integrated
# VIS and NIR
# total
# snowpack
# absorption:
# abs_vis = sum(flx_slr(1:vis_max_idx).*(1 - albedo(1:vis_max_idx)));
# abs_nir = sum(flx_slr(vis_max_idx + 1:nir_max_idx).*(1 - albedo(vis_max_idx + 1:nir_max_idx)));

# % Output
# diagnostics:
# % 1.
# The
# basics:
# data_out(:, 1) = wvl; % spectral
# wavelength
# bands(um)
# data_out(:, 2) = albedo; % spectral
# hemispheric
# albedo
# data_out(1, 3) = alb; % solar
# broadband
# albedo
# data_out(2, 3) = alb_vis; % visible(0.3 - 0.7
# um) albedo
# data_out(3, 3) = alb_nir; % near - IR(0.7 - 5.0
# um) albedo
# data_out(4, 3) = sum(F_abs_slr); % total
# radiative
# absorption
# by
# all
# snow
# layers(not including
# underlying
# substrate)
# data_out(5, 3) = F_abs_slr(1); % top
# layer
# solar
# absorption
# if (nbr_lyr > 1)
# data_out(6, 3) = F_abs_slr(2); % 2
# nd
# layer
# solar
# absorption
# else
# data_out(6, 3) = NaN; % 2
# nd
# layer
# solar
# absorption
# end;
# data_out(:, 5)  = sum(F_abs, 2); % spectral
# absorption

# % more
# detail:
# if (BND_TYP == 1)
#         % different band-weighted albedos:
#     data_out(1, 4) = alb1;
#     data_out(2, 4) = alb2;
#     data_out(3, 4) = alb3;
#     data_out(4, 4) = alb4;
#     data_out(5, 4) = alb5;
#     data_out(7, 3) = alb6;
#     data_out(8, 3) = alb7;
#     data_out(9, 3) = alb8;
#     end;

#     data_out(6, 4) = F_abs_slr(1); %
# top
# layer
# solar
# absorption
# data_out(7, 4) = F_abs_vis(1); % top
# layer
# VIS
# absorption
# data_out(8, 4) = F_abs_nir(1); % top
# layer
# NIR
# absorption

# if (nbr_lyr > 1)
# data_out(9, 4) = F_abs_slr(2);
# data_out(10, 4) = F_abs_vis(2);
# data_out(11, 4) = F_abs_nir(2);
# if (nbr_lyr > 2)
# data_out(12, 4) = F_abs_slr(3);
# data_out(13, 4) = F_abs_vis(3);
# data_out(14, 4) = F_abs_nir(3);
# if (nbr_lyr > 3)
# data_out(15, 4) = F_abs_slr(4);
# data_out(16, 4) = F_abs_vis(4);
# data_out(17, 4) = F_abs_nir(4);
# if (nbr_lyr > 4)
# data_out(18, 4) = F_abs_slr(5);
# data_out(19, 4) = F_abs_vis(5);
# data_out(20, 4) = F_abs_nir(5);
# end;
# end;
# end;
# end;

# data_out(18, 4) = F_abs_btm; % solar absorption by underlying surface
# data_out(19, 4) = F_abs_vis_btm; % VIS absorption by underlying surface
# data_out(20, 4) = F_abs_nir_btm; % NIR absorption by underlying surface
# data_out(21, 4) = sum((mu_not * pi * Fs))+sum(Fd); % total downwelling energy on upper boundary
# data_out([1:5], 6) = heat_rt; % JC
# EDIT
# output
# radiative
# heating
# rate in each
# layer in K / hr
# data_out([1: 470], [7: 11]) = intensity2; % JC
# EDIT
# output
# planar
# intensity
# per
# layer

# % plot
# modeled
# albedo:
# if (1 == 0)
# plot(wvl, albedo, 'k-');
# axis([0.3 2.5 0 1]);
# grid on;
# end;

