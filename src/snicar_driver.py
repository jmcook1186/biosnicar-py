#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("./src")
from setup_snicar import *
from classes import *
from column_OPs import *
from biooptical_funcs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver
from validate_inputs import *
from display import *

###################
# BIO-OPTICAL MODEL
###################

# optionally run the bio-optical modelk to add new impurity optical properties to
# the BioSNICAR database. Commewnted out by default as we expect our default lap
# database to be sufficient for most users.

# run_biooptical_model()

###########################
# RADIATIVE TRANSFER MODEL
###########################

# define input file
input_file = "./src/inputs2.yaml"

# first build classes from config file and validate their contents
(
    ice,
    illumination,
    rt_config,
    model_config,
    plot_config,
    impurities,
) = setup_snicar(input_file)

# validate inputs to ensure no invalid combinations have been chosen
status = validate_inputs(ice, rt_config, model_config, illumination, impurities)

# now get the optical properties of the ice column
ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
tau, ssa, g, L_snw = mix_in_impurities(
    ssa_snw, g_snw, mac_snw, ice, impurities, model_config
)

# now run one or both of the radiative transfer solvers
outputs1 = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config)

#outputs2 = toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)

# plot and print output data
plot_albedo(plot_config, model_config, outputs1.albedo)
display_out_data(outputs1)
