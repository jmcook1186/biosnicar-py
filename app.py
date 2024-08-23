import numpy as np
import pandas as pd
import biosnicar
from biosnicar.adding_doubling_solver import adding_doubling_solver
from biosnicar.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.setup_snicar import setup_snicar
from biosnicar.validate_inputs import validate_inputs
import plotly.express as px
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
snow ❄️ ice 🧊 and life :space_invader: albedo model (v{biosnicar.__version__})

[GitHub](https://github.com/jmcook1186/biosnicar-py)
[Documentation](https://biosnicar.vercel.app)

*Note that impurities are assumed to exist in the upper 2 cm of the snow or ice only, and that the grain shape in the granular layer set up is spherical.
To access other configurations, download and run the full model as Python code instead.*
"""
)
st.markdown("""---""")


st.sidebar.header("Ice")
layer = st.sidebar.selectbox("Layer Type", ("grains", "solid ice"))
thickness = st.sidebar.number_input(
    "Thickness (meters)", 0.0, 10.0, 3.0, help="ice layer thickness in meters"
)
radius = st.sidebar.number_input(
    "Radius (microns)",
    0,
    10000,
    5000,
    1000,
    help="grain radius at layer type grain or bubble radius Layer type solid in each layer",
)

density = st.sidebar.number_input("Density (kg/m^3)", 0, 1000, 700)
lwc = st.sidebar.number_input("Liquid water content", 0.0, 1.0, 0.01)

st.sidebar.header("Impurities")
black_carbon = st.sidebar.number_input("Black carbon conc (ppb)", 0, 1000000, 0)
glacier_algae = st.sidebar.number_input("Glacier algae conc (cells/mL)", 0, 1000000, 0)
snow_algae = st.sidebar.number_input("Snow algae conc (cells/mL)", 0, 1000000, 0)

st.sidebar.header("Sun")
solar_zenith_angle = st.sidebar.number_input("Solar zenith angle (degrees)", 1, 89, 50)


def run_snicar(
    layer: str,
    thickness: float,
    radius: int,
    density: int,
    lwc: float,
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

    input_file = "app_inputs.yaml"

    if layer == "grains":
        layer = 4
    elif layer == "solid ice":
        layer = 1
    else:
        raise ValueError("invalid layer type")

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
    ice.dz = [0.02, thickness-0.02]
    ice.rds = [radius, radius]
    ice.rho = [density, density]
    ice.lwc = [lwc, lwc]

    impurities[0].conc = [black_carbon, 0]
    impurities[1].conc = [snow_algae, 0]
    impurities[2].conc = [glacier_algae, 0]

    illumination.solzen = solar_zenith_angle
    illumination.calculate_irradiance()
    # validate inputs to ensure no invalid combinations have been chosen
    status = validate_inputs(ice, illumination, impurities)

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
    albedo.index.name = "wavelength (microns)"
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
        labels={"index": "wavelengths (microns)", "value": "Albedo"},
    )
    fig.update_layout(showlegend=False)
    return fig


result = run_snicar(
    layer,
    thickness,
    radius,
    density,
    lwc,
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
