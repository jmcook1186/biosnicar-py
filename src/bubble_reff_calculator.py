"""
This script is for taking specific surface area for bubbly ice an calculating the
effective radius of the bubbles assuming a lognormal bubble size distribution.

The reason this is required is that the well-known conversion between SSA and r_eff:

r_eff = 3/(P_i * SSA) where P_i is density of ice (917 kg m-3)

gives the effective radius of a discrete grain of given SSA, for a collection of bubbles
in a bulk medium of ice a more nuanced calculation is required.

The derivation of the calculations are explained in the pdf ./Assets/SSA_derivation.pdf
from Chloe Whicker (UMich).

For initial experiments in simulating GrIS glacier ice beneath the weathered crust
I have taken SSA values measured in ice cores of bubbly glacier ice from the
Allan Hills (Antarctica)
by CT-scanning by Dadic et al (2013), in the absence of data from Greenland.
"""

import numpy as np


# SET CONSTANTS
RHO_ICE = 917  # density of pure ice
RHO_AIR = 1.025  # density of air
SIGMA_G = 1.5  # unitless
SIGMA_TILDE_G = np.log(SIGMA_G)


# SET VARIABLES
# Dadic et al. densities from calliper measurements
# rho = [460, 668, 777, 864, 894] # kg / m3
# Dadic et al. densities from Micro-CT
rho = [460, 737, 820, 891, 894]  # kg / m3
# Dadic et al. speciic surface area from Micro-CT
# ssa = [15.4, 3.8, 1.9, 0.39, 0.16]
# Dadic et al. specific surface area from RT model
ssa = [15.4, 7.8, 3.9, 0.65, 0.16]  # Specific Surface Area units: m2/kg


# SETUP EMPTY ARRAYS
D_eff = np.zeros(len(rho))

# CALCULATE D_eff IN LOOP
for i in np.arange(0, len(rho), 1):

    if rho[i] < 550:  # anything with a density smaller than 550 is represented as snow

        D_eff[i] = 6 / (RHO_ICE * ssa[i])  # D_eff in meters

    else:
        V_air = (rho[i] - RHO_ICE) / (-RHO_ICE + RHO_AIR)  # unitless volume fraction
        D_eff[i] = (6 * V_air) / (rho[i] * ssa[i])  # D_eff in meters

# Number of bubbles
No = (6 * V_air) / (np.pi * D_eff**3) * np.exp(3 * SIGMA_TILDE_G**2)

# convert effective diameter in meters to effective radius in microns
r_eff_micron = D_eff / 2 * 1e6

print(np.round(r_eff_micron, 0))
