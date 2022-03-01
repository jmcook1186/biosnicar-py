#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
(
    ice,
    illumination,
    rt_config,
    model_config,
    plot_config,
    impurities,
) = setup_snicar()

status = validate_inputs(ice, rt_config, model_config, illumination, impurities)


# now get the optical properties of the ice column
ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
tau, ssa, g, L_snw = mix_in_impurities(
    ssa_snw, g_snw, mac_snw, ice, impurities, model_config
)

# now run one or both of the radiative transfer solvers
outputs1 = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config)

outputs2 = toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)
# plot_albedo(plot_config,model_config, outputs2.albedo)
print(outputs1.BBA)

plot_albedo(plot_config, model_config, outputs1.albedo)
display_out_data(outputs1)
