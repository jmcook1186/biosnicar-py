import numpy as np
import pandas as pd
import plotly.express as px
import biosnicar
from biosnicar.adding_doubling_solver import adding_doubling_solver
from biosnicar.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.setup_snicar import setup_snicar
from biosnicar.validate_inputs import validate_inputs

import streamlit as st

st.set_page_config(
    page_title="biosnicar",
    page_icon="❄️",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://biosnicar-go-py.readthedocs.io/en/latest/",
        "Report a bug": "https://github.com/jmcook1186/biosnicar-py/issues",
        "About": f"# biosnicar frontent version {biosnicar.__version__}",
    },
)

st.title(f"biosnicar")
st.markdown(
    f"""
snow ❄️ and ice 🧊 albedo model (v{biosnicar.__version__})

[GitHub](https://github.com/jmcook1186/biosnicar-py)
[Documentation](https://biosnicar-go-py.readthedocs.io/en/latest/)
"""
)
st.markdown("""---""")


st.sidebar.header("Ice")
layer = st.sidebar.selectbox("Layer Type", ("grains", "solid"))
thickness = st.sidebar.number_input(
    "Thickness (meters)", 0.0, 10.0, 3.0, help="ice layer thickness in meters"
)
radius = st.sidebar.number_input(
    "Radius (microns)",
    0,
    10000,
    7000,
    help="grain radius at layer type grain or bubble radius Layer type solid in each layer",
)
density = st.sidebar.number_input("Density (kg/m^3)", 0, 1000, 700)

st.sidebar.header("Impurities")
black_carbon = st.sidebar.number_input("Black carbon conc (ppb)", 0, 10000, 0)
glacier_algae = st.sidebar.number_input("Glacier algae conc (cells/mL)", 0, 10000, 0)
snow_algae = st.sidebar.number_input("Snow algae conc (cells/mL)", 0, 10000, 0)

st.sidebar.header("Sun")
solar_zenith_angle = st.sidebar.number_input("Solar zenith angel (degree)", 1, 89, 50)


def run_snicar(
    layer: str,
    thickness: float,
    radius: int,
    density: int,
    black_carbon: int,
    glacier_algae: int,
    snow_algae: int,
    solar_zenith_angle: int,
) -> dict:
    """Runs biosnicar

    Args:
        layer (str): _description_
        thickness (float): _description_
        radius (int): _description_
        density (int): _description_
        black_carbon (int): _description_
        glacier_algae (int): _description_
        snow_algae (int): _description_
        solar_zenith_angle (int): _description_

    Returns:
        dict: Dict with result for display.
    """

    input_file = "app/api/inputs.yaml"

    if layer == "grains":
        layer = 0
    elif layer == "solid":
        layer = 1

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
    rounded_broandband_albedo = np.round(outputs.BBA, 2)
    wave_length = np.arange(0.205, 5, 0.01)
    albedo = pd.Series(outputs.albedo, index=wave_length, name="albedo")
    albedo.index.name = "wave lenght (microns)"
    albedo_csv = albedo.to_csv()

    return {
        "albedo": albedo,
        "broadband": rounded_broandband_albedo,
        "status": status,
        "albedo_csv": albedo_csv,
    }


def plot_albedo(albedo: pd.Series):
    fig = px.line(
        result["albedo"],
        range_y=[0, 1],
        range_x=[0.205, 2.5],
        labels={"index": "wave lenghts (microns)", "value": "Albedo"},
    )
    fig.update_layout(showlegend=False)
    return fig


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

# display results
st.metric("Broadband Albedo", result["broadband"])
st.plotly_chart(plot_albedo(result["albedo"]))
st.download_button("download data", data=result["albedo_csv"], file_name="albedo.csv")

with st.expander("Show raw data"):
    st.dataframe(result["albedo"])
