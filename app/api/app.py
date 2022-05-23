#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script configures and deploys the BioSNICAR web app.

The BioSNICAR web-app offers convenient browser access to the
BioSNICAR model code via a simple GUI. It has reduced functionality
compared to the full BioSNICAR model but is likely to be sufficient
for many users.

Here are two functions, root() and run_snicar(). The former sets up
a default route for the frontend to connect to. It renders index.html
to the browser, where the user can define their input values.
These are passed as a POST request to run_snicar() which
saves albedo data to file (jpg and csv). The frontend detects the
change in file contents and re-renders the browser contents.

"""
import time

import matplotlib.pyplot as plt
import numpy as np
import requests
from biosnicar.adding_doubling_solver import adding_doubling_solver
from biosnicar.biooptical_funcs import *
from biosnicar.classes import *
from biosnicar.column_OPs import *
from biosnicar.display import *
from biosnicar.setup_snicar import *
from biosnicar.toon_rt_solver import toon_solver
from biosnicar.validate_inputs import *
from flask import Flask, Response, request
from flask_cors import CORS, cross_origin

app = Flask(__name__, static_folder="../build", static_url_path="/")

# set root url for host as a dynamic variable
port = int(os.environ.get("PORT", 5000))
host = "http://localhost:5000/"
success = False

# enable cross domain requests
# because flask is on port 5000, react on port 3000
app = Flask(__name__)
cors = CORS(app)
app.config["CORS_HEADERS"] = ["Content-Type"]


@app.route("/")
def root():
    """
    Sets up default route to connect backend to frontend.

    Args:
        None

    Returns:
        renders index.html to browser
    """
    return app.send_static_file("index.html")


# the actual func to run
@app.route("/app/model", methods=["POST"])
@cross_origin()
def run_snicar():
    """
    Runs BioSNICAR model and renders output to webpage.

    This function takes the use inputs from the POST request submitted
    via the React frontend. These inputs are passed to the BioSNCAR
    model and a plot of spectral albedo is saved to file.
    The broadband albedo, rounded to 2 decimal places, is printed on the
    figure.

    Args:
        None, but receives args as POSt http request via
        /run_snicar route.

    Returns:
        res: http response code

    """

    input_file = "app/api/inputs.yaml"

    lyr_typ = int(request.json["lyr_typ"])
    dz = float(request.json["dz"])
    r_eff = int(request.json["r_eff"])
    rho = int(request.json["rho"])
    bc = int(request.json["bc"])
    glacier_algae = int(request.json["glacier_algae"])
    snow_algae = int(request.json["snow_algae"])
    zenith = int(request.json["zenith"])

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

    impurities[0].conc = [bc, 0]
    impurities[1].conc = [snow_algae, 0]
    impurities[2].conc = [glacier_algae, 0]

    illumination.solzen = zenith
    illumination.calculate_irradiance()
    # validate inputs to ensure no invalid combinations have been chosen
    status = validate_inputs(ice, rt_config, model_config, illumination, impurities)

    # now get the optical properties of the ice column
    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(
        ssa_snw, g_snw, mac_snw, ice, impurities, model_config
    )

    # now run one or both of the radiative transfer solvers
    outputs = adding_doubling_solver(
        tau, ssa, g, L_snw, ice, illumination, model_config
    )
    rounded_BBA = np.round(outputs.BBA, 2)
    wvl = np.arange(0.205, 5, 0.01)
    plt.figure()
    plt.plot(wvl, outputs.albedo, linewidth=1)
    plt.ylim(0, 1)
    plt.xlim(0.2, 2.5)
    plt.text(1.5, 0.8, f"broadband albedo\n{rounded_BBA}", fontsize=20)
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Albedo")
    plt.grid(False)

    plt.tight_layout()
    plt.savefig("app/src/outputs/albedo.jpg")
    plt.close()

    np.savetxt("app/src/outputs/albedo.csv", outputs.albedo, delimiter=",")

    res = Response("success")
    res.headers["Access-Control-Allow-Origin"] = "*"

    return res


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=port, debug=True)
