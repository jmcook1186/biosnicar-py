"""
BioSNICAR_GO DRIVER SCRIPT

This script is used to run the 2-stream radiative transfer
model BioSNICAR_GO. Here solvers are called and the results plotted.

NB: if using only granular layers, recommend using the faster toon et al
tridiagonal matix solver (by setting toon = True), however this will not
include any specular reflection components. If solid ice layers are
included in the ice/snow column, the ADDING-DOUBLING solver must be used
(i.e. add_double = True).

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
from validate_inputs import *
from display import *


# first build classes from config file and validate their contents
ice, illumination, rt_config, model_config, plot_config, display_config, impurities = setup_snicar()
status = validate_inputs(ice, rt_config, model_config, illumination, impurities)

# now get the optical properties of the ice column
ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
tau, ssa, g, L_snw = mix_in_impurities(
    ssa_snw, g_snw, mac_snw, ice, impurities, model_config
)

# now run one or both of the radiative transfer solvers
outputs = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config)

# outputs=toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)
# plot_albedo(plot_config,model_config, outputs2.albedo)
print(outputs.BBA)

plot_albedo(plot_config, model_config, outputs.albedo)
display_out_data(outputs)