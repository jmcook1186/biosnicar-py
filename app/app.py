
import sys
sys.path.append("./src")
from flask import Flask, request, render_template, redirect, url_for
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
app = Flask(__name__, static_url_path = "/templates", static_folder = "templates")

# set root url for host as a dynamic variable
host = 'http://localhost:5000/'
success = False
####################   ROUTE 1: /welcome   #############################
#### renders welcome page, sends POST request to /classifier route #####
#####    and redirects to /output when POST request completed      #####

@app.route('/welcome', methods=['POST', 'GET'])
def welcome():
    # this function creates a homepage for the app - when the browser is
    # directed to localhost:5000/welcome it displays a welcome message and
    # an example image stored as a static file in the /templates folder

    global success

    # react to data submitted from web form in homepage.html
    if request.method == "POST":
        # get data submitted from html form
        lyr_typ = request.form['lyr_typ']
        dz = request.form['dz']
        r_eff = request.form['r_eff']
        rho = request.form['rho']


        # create POST request using requests package
        api_url = str(host + 'model') # set url to send request to
        payload = {'lyr_typ': lyr_typ, 'dz': dz, 'r_eff': r_eff, 'rho': rho}
        
        
        # send request with payload = create_row_data
        r = requests.post(url=api_url, json=payload)
        print("PRINTING REQUEST STATUS: CODE = {}, REASON = {}, TEXT = {} ".format(r.status_code, r.reason, r.text))

        if success == True:
            return redirect(url_for('output')) # on request success redirect browser to URL for output() function
        # else:
        #     return redirect(url_for('output_fail'))

    return render_template('welcome.html', title='BioSNICAR') # while post request incomplete render homepage


######################   ROUTE 2: /output   ###########################
#######     renders responsepage html template via /output     ########

@app.route("/output")
def output():
    return render_template("response.html", title="BioSNICAR")


# @app.route("/output_fail")
# def output_fail():
#     return render_template("responsepage_fail.html", title="Ice Classifiers")


# @app.route("/info_page")
# def info_page():
#     return render_template("info_page.html", title="Ice Classifiers")

####################   ROUTE 3: /model   #############################
#### this route contains the model calls. Input data is gathered from  ####
#### the html form in homepage.html & parsed to the successive python  ####
#### scripts comprising the IceClassifiers model. Output image saves   ####
#### to ./templates to make the file available as a static resource to ####
#### display in the browser when responsepage.html renders via /output ####

@app.route('/model', methods=['POST'])
def post_function():

    global success
    json = request.get_json()
    lyr_typ = int(json["lyr_typ"])
    dz = float(json["dz"])
    r_eff = int(json["r_eff"])
    rho = int(json["rho"])

    input_file = "app/inputs.yaml"

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
    wvl = np.arange(0.205, 5, 0.01)
    plt.plot(wvl, outputs.albedo)
    plt.ylim(0,1)
    plt.xlim(0.2,2.5)
    plt.savefig("app/templates/outputs/albedo.jpg")


    success=True

    return


# run app
if __name__ == "__main__":
    app.run(debug=True)






