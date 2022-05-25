import streamlit as st
import numpy as np
import pandas as pd
from biosnicar.setup_snicar import setup_snicar
from biosnicar.validate_inputs import validate_inputs
from biosnicar.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.adding_doubling_solver import adding_doubling_solver
import matplotlib.pyplot as plt
import plotly.express as px


st.title("biosnicar")
st.write("snow and ice albedo model")

layer = st.sidebar.selectbox("Layer Type", (0, 1))
thickness = st.sidebar.number_input("Thickness (meters)", 0.0, 10.0, 3.0)
radius = st.sidebar.number_input("Radius (microns)", 0, 10000, 7000)
density = st.sidebar.number_input("Density (kg/m^3)", 0, 1000, 700)
black_carbon = st.sidebar.number_input("Black carbon conc (ppb)", 0, 1000, 0)
glacier_algae = st.sidebar.number_input("Glacier algae conc (cells/mL)", 0, 1000, 0)
snow_algae = st.sidebar.number_input("Snow algae conc (cells/mL)", 0, 1000, 0)
solar_zenith_angle = st.sidebar.number_input("Solar zenith angel (degree)", 1, 89, 60)


def run_snicar(
    layer: int,
    thickness: float,
    radius: int,
    density: int,
    black_carbon: int,
    glacier_algae: int,
    snow_algae: int,
    solar_zenith_angle: int,
):
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
    ice.layer_type = [layer, layer]
    ice.dz = [0.001, thickness]
    ice.rds = [radius, radius]
    ice.rho = [density, density]

    impurities[0].conc = [black_carbon, 0]
    impurities[1].conc = [snow_algae, 0]
    impurities[2].conc = [glacier_algae, 0]

    illumination.solzen = solar_zenith_angle
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
    wave_length = np.arange(0.205, 5, 0.01)
    albedo = pd.Series(outputs.albedo, index=wave_length, name="albedo")

    return {"albedo": albedo, "broadband": rounded_BBA}


st.write("You selected:", layer)
result = run_snicar(
    layer,
    thickness,
    radius,
    density,
    black_carbon,
    glacier_algae,
    snow_algae,
    solar_zenith_angle,
)
st.metric("Broadband Albedo", result["broadband"])
st.plotly_chart(
    px.line(
        result["albedo"], range_y=[0, 1], labels=["Wavelenghts (microns)", "Albedo"]
    )
)
st.dataframe(result["albedo"])
