import numpy as np
import pandas as pd
import biosnicar
from biosnicar.drivers.run_model import run_model
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
    """Runs biosnicar via the unified run_model() entry point.

    Args:
        layer: "grains" or "solid ice"
        thickness: total ice thickness in metres
        radius: grain/bubble radius in microns
        density: layer density in kg/m³
        lwc: liquid water content
        black_carbon: BC concentration in ppb
        glacier_algae: glacier algae concentration in cells/mL
        snow_algae: snow algae concentration in cells/mL
        solar_zenith_angle: solar zenith angle in degrees

    Returns:
        dict with albedo Series, broadband albedo, and CSV data.
    """
    if layer == "grains":
        layer_type = 3
    elif layer == "solid ice":
        layer_type = 1
    else:
        raise ValueError("invalid layer type")

    outputs = run_model(
        solver="adding-doubling",
        validate=True,
        solzen=solar_zenith_angle,
        rds=radius,
        rho=density,
        dz=[0.02, thickness - 0.02],
        lwc=lwc,
        layer_type=layer_type,
        impurity_0_conc=[black_carbon, 0],
        impurity_1_conc=[snow_algae, 0],
        impurity_2_conc=[glacier_algae, 0],
    )

    rounded_broandband_albedo = np.round(outputs.BBA, 2)
    wave_length = np.arange(0.205, 5, 0.01)
    albedo = pd.Series(outputs.albedo, index=wave_length, name="albedo")
    albedo.index.name = "wavelength (microns)"
    albedo_csv = albedo.to_csv()

    return {
        "albedo": albedo,
        "broadband": rounded_broandband_albedo,
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
