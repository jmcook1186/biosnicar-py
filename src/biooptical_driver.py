#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joe, lou

Driver for the bio-optical model developed by Cook et al. 2017, 2020 to
calculate algae optical properties and save them to a netcdf files directly
usable in BioSNICAR. 

Note:



"""

from biooptical_funcs import *
from classes import *
import sys
sys.path.append("./src")

# --------------------------------------------------------------------------------------
# CALCULATIONS OF ABSORPTION PROPERTIES
# --------------------------------------------------------------------------------------

bio_optical_config = BioOpticalConfig()
ACS = get_absorption_cross_section(bio_optical_config)
k = calculate_k(bio_optical_config, ACS)
wvl_rescaled_BioSNICAR, ACS_rescaled_BioSNICAR,\
k_rescaled_BioSNICAR, n_rescaled_BioSNICAR = rescale_480band(bio_optical_config, 
                                                              ACS, k)
plot_k_n_ACS(bio_optical_config, ACS, k)

# # --------------------------------------------------------------------------------------
# # CALCULATIONS OF SCATTERING PROPERTIES
# # --------------------------------------------------------------------------------------

assym, ss_alb = calculate_ssps(bio_optical_config, k_rescaled_BioSNICAR, 
                                wvl_rescaled_BioSNICAR,n_rescaled_BioSNICAR)


# # --------------------------------------------------------------------------------------
# # SAVING DATA IN NETCDF
# # --------------------------------------------------------------------------------------

net_cdf_updater(bio_optical_config, assym, ss_alb, ACS_rescaled_BioSNICAR, wvl_rescaled_BioSNICAR)


