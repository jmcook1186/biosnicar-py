#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script configures and deploys the BioSNICAR web app. 

The BioSNICAR web-app offers convenient browser access to the
BioSNICAR model code via a simple GUI. It has reduced functionality
compared to the full BioSNICAR model but is likely to be sufficient
for many users.

Here are two functions, main() and model(). The former sets up
a landing page with a simpel Flask form for parsing user input
values. These are passed as a POST request to model() which
saves albedo data to file. This file is rendered in a response
page that is redirected to automatically.

"""
import sys
sys.path.append("../../src")
sys.path.append("./src")
sys.path.append("../app")
from flask import Flask, request, Response
from flask_cors import CORS, cross_origin
import time
import numpy as np
import requests
import matplotlib.pyplot as plt
from setup_snicar import *
from classes import *
from column_OPs import *
from biooptical_funcs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver
from validate_inputs import *
from display import *


# Flask static paths set to /templates to enable retrieval of default image by index.html
app = Flask(__name__, static_folder="/")

# set root url for host as a dynamic variable
port = int(os.environ.get("PORT", 5000))
host = 'http://localhost:5000/'
success = False

# enable cross domain requests
# because flask is on port 5000, react on port 3000
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = ['Content-Type']


# the actual func to run
@app.route('/model', methods=['POST'])
@cross_origin()
def run_snicar():
    """
    Runs BioSNICAR model and renders output to webpage.

    This function takes the use inputs from the POST request submitted
    by the user via the webform on route "/". These inputs are passed
    to the BioSNCAR model and a plot of spectral albedo is saved to file.
    The broadband albedo, rounded to 2 decimal places, is printed on the
    figure. A successful model run triggers an auto redirect to 
    the /output page where the figure is rendered.

    """

    input_file = "react-app/app/inputs.yaml"

    lyr_typ = int(request.json['lyr_typ'])
    dz = float(request.json['dz'])
    r_eff = int(request.json['r_eff'])
    rho = int(request.json['rho'])

    print(lyr_typ)
    print(dz)
    print(r_eff)
    print(rho)

    # first build classes from config file and validate their contents
    (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
        impurities,
    ) = setup_snicar(input_file)


    # load base classes from default inputs.yaml
    # then adjust for user inputs
    ice.layer_type = [lyr_typ, lyr_typ]
    ice.dz = [0.001, dz]
    ice.rds = [r_eff, r_eff]
    ice.rho = [rho, rho]


    # validate inputs to ensure no invalid combinations have been chosen
    status = validate_inputs(ice, rt_config, model_config, illumination, impurities)

    # now get the optical properties of the ice column
    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(
        ssa_snw, g_snw, mac_snw, ice, impurities, model_config
    )

    # now run one or both of the radiative transfer solvers
    outputs = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config)
    rounded_BBA = np.round(outputs.BBA, 2)
    wvl = np.arange(0.205, 5, 0.01)
    plt.plot(wvl, outputs.albedo)
    plt.ylim(0,1)
    plt.xlim(0.2,2.5)
    plt.text(1.75, 0.85, f"BBA: {rounded_BBA}", fontsize = 20)
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Albedo")
    plt.tight_layout()
    plt.savefig("react-app/src/outputs/albedo.jpg")
    plt.close()

    res = Response("success")
    res.headers['Access-Control-Allow-Origin'] = '*'

    return res


if __name__ == '__main__':
    app.run(host = '0.0.0.0', port = port, debug=True)






