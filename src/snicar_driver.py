"""
BioSNICAR_GO DRIVER SCRIPT

This script is used to configure the 2-stream radiative transfer
model BioSNICAR_GO. Here variable values are defined, the model called
and the results plotted.

NB. Setting Mie = 1, GO = 0 and algal impurities = 0 is equivalent to
running the original SNICAR model of Flanner et al. (2007, 2009)

NB: if using only granular layers, recommend using the faster toon et al
tridiagonal matix solver (by setting toon = True), however this will not
include any specular reflection components. If solid ice layers are
included in the ice/snow column, the ADDING-DOUBLING solver must be used
(i.e. add_double = True).

glacier algae MAC files are provided in units of m2/cell, which requires
a unit conversion that is not applied to the other LAPS. The conversion
is selectivlely applid by indexing to the last element in the LAP
lists, meaning glacier algae must always be the final LAP, even if more
LAPs are added in future.

Author: Joseph Cook, June 2021
"""

import collections as c
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from setup_snicar import *
from classes import *
from column_OPs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver


ice = Ice()
illumination = Illumination()
rt_config = RTConfig()
model_config = ModelConfig()

impurities = build_impurities_array()

for i in [0,  50000]:
    
    impurities[0].conc = [i,0]

    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(ssa_snw, g_snw, mac_snw, ice, impurities, model_config)
    
    outputs1 = toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)
    outputs2 = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)

    # plt.figure()
    # plt.plot(outputs1.albedo)
    # plt.show()