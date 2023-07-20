#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pathlib import Path
from validate_inputs import validate_inputs
from adding_doubling_solver import adding_doubling_solver
from column_OPs import get_layer_OPs, mix_in_impurities
from display import display_out_data, plot_albedo
from setup_snicar import setup_snicar
from toon_rt_solver import toon_solver


def get_albedo(solver, plot, validate):
    (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
        impurities,
    ) = setup_snicar()

    if validate:
        validate_inputs(ice, illumination, impurities)

    # now get the optical properties of the ice column
    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(
        ssa_snw, g_snw, mac_snw, ice, impurities, model_config
    )
    # now run one or both of the radiative transfer solvers
    if solver == "toon":
        print("\nRunning biosnicar with the Toon solver\n")
        outputs = toon_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )
    elif solver == "adding-doubling":
        print("\nRunning biosnicar with the adding-doubling solver\n")
        outputs = adding_doubling_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config
        )
    else:
        return "solver not recognized, please choose toon or adding-doubling"

    if plot:
        plot_albedo(plot_config, model_config, outputs.albedo)
    display_out_data(outputs)
    return outputs.albedo
