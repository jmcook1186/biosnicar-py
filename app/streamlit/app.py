import streamlit as st
import biosnicar


st.title("biosnicar")
st.write("snow and ice albedo model")

layer = st.sidebar.selectbox("Layer Type", (0, 1, 2))
thickness = st.sidebar.number_input("Thickness (meters)", 0.0, 10.0)
radius = st.sidebar.number_input("Radius (microns)", 0.0, 10000.0)
density = st.sidebar.number_input("Density (kg/m^3)", 0.0, 1000.0, 700.0)
black_carbon = st.sidebar.number_input("Black carbon conc (ppb)", 0.0, 1000.0, 0.0)
glacier_algae = st.sidebar.number_input(
    "Glacier algae conc (cells/mL)", 0.0, 1000.0, 0.0
)
snow_algae = st.sidebar.number_input("Snow algae conc (cells/mL)", 0.0, 1000.0, 0.0)
solar_zenit_angle = st.sidebar.number_input(
    "Solar zenith angel (degree)", 0.0, 90.0, 0.0
)

st.write("You selected:", layer)
