#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from biosnicar import run_model

outputs = run_model(
    solver="adding-doubling",
    plot=False,
    validate=True,
    solzen=50,
    layer_type=[1, 1, 1, 1, 1],
    dz=[0.02, 0.05, 0.05, 0.05, 0.83],
    rds=[800, 900, 1000, 1100, 1200],
    rho=[500, 600, 700, 750, 800],
    lwc=[0, 0, 0, 0, 0],
    impurity_0_conc=[0, 0, 0, 0, 0],       # black carbon (ppb)
    impurity_1_conc=[40000, 10000, 0, 0, 0],  # snow algae (cells/mL)
    impurity_2_conc=[0, 0, 0, 0, 0],       # glacier algae (cells/mL)
)

print(f"Broadband albedo: {outputs.BBA}")
