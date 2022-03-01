# BioSNICAR

<img src="./assets/example-output.jpg" width=500>


## Introduction

BioSNICAR is a set of Python scripts that predict the spectral albedo of snow and glacier ice between 200nm to 5000nm given information about the illumination conditions, ice structure and the type and concentration of light absorbing particulates (LAPs) externally mixed with the snow/ice. The jumping off point for this model was legacy FORTRAN and Matlab code of the SNICAR model developed by Flanner et al. (2007), which solves the radiative transfer equations simplified into a two-stream scheme after Toon et al. 1989. Two solvers are available: the original SNICAR matrix solver typically representing ice and snow with grains (Toon et al. 1989) and the Adding-Doubling (AD) solver representing the ice as a solid medium with air bubbles and allowing the incorporation of Fresnel reflecting layers (Brieglib and Light, 2007, Dang et al. 2019, Wicker et al. 2022). BioSNICAR couples SNICAR to a bio-optical model that allows for calculation of optical properties of snow and glacier algae to load into the model as LAPs (Cook et al. 2017, 2020). This functionality, along with the vectorized AD solver formulation, accessible user interface and applicability to a very wide range of surface conditions are the unique selling points of this implementation. This code is also very actively maintained and we welcome contributions from the community to help make BioSNICAR a useful tool for a diverse range of cryosphere scientists.

## How to use

### Installing Environment/Dependencies

If you do not have Python installed, download Python >3.6. It is recommended to use a fresh environment using conda or venv. Once activated, install the project dependencies with:

```
pip install -r requirements.txt

```

Finally, if you do not wish to install anything on your computer, but you use VSCode and Docker, then you can use the devcontainer config provided to run this code in a remote container. This requires the "remote containers" extension to be added to VSCode. Further instructions are available here: https://code.visualstudio.com/docs/remote/containers


### Running the model

The model driver and all the core source code can be found in `/src`. From the top level directory (`~/BioSNICAR_GO_PY`), run:

`python ./src/snicar_driver.py`

This will run the model with all the default settings. The user will see a list of output values printed to the console and a spectral albedo plot appear in a separate window. The code can also be run in an interactive session (Jupyter/iPython) in which case the relevant data and figure will appear in the interactive console. 

Most users will want to experiment with changing input parameters. This is achieved by adjusting the values in the config file `inputs.yaml`. The nature of each parameter is described in in-line annotations to guide the user. Invalid combinations of values will be rejected by our error-checking code. Most users should have no reason to modify any other file in this repository except for the those in `inputs.yaml`.

More complex applications of the model code, for example model inversions, field/model comparisons etc are included under `/experiments`, with details provided in that module's own README.

We have also maintained a separate version of the BioSNICAR codebase that uses a "functional" prorgamming style rather than the object-oriented approach taken here. We refer to this as BioSNICAR Classic and it is available in the `classic` branch of this repository. it might be useful for people already familiar with FORTRAN or Matblab implementations from previous literature. The two branches are entirely equivalent int heir simulations but very different in their programmign style. The object oriented approach is preferred because it is more Pythonic, more flexible and easier to debug.

## Choosing Inputs
It is straightforward to adjust the model configuration by updating the values in `inputs.yaml`. However there is a lot of nuance to setting up the model to provide realistic simulations, and the meaning of the various parameters is not always obvious. for this reason we have put together a guide. Please refer to `inputs.md` for explanations of each parameter. 

## Theoretical Background

A detailed description of the theory on which the model is based can be found in Flanner et al. 2021. Briefly, the spectral albedo is computed from the single scattering  properties of snow/ice and the light absorbing particulates (LAPs), along with the concentrations of LAPs and the physical properties of snow/ice (density, grain/pore size). These can vary in each layer, allowing for representation of e.g. porous weathered crust above unweathered ice with an uppermost mm of liquid water containing algae, or deep layers of ancient dust. Below are detailed the different options to calculate the single scattering properties.

### Snow and ice 

The snow and ice single scattering properties are generated from the refractive index of ice from Warren 1984, Warren 2008 or Picard 2016 and all parametrisation is done in the driver of the model. Snow and ice are represented either as a collection of grains, or as a continuous layer with embedded air bubbles. This second option requires to use the AD solver with fresnel layers (ADD_DOUBLE and LAYER_TYPE toggled on in the model driver). The embedded bubbles are represented as spheres of air surrounded by ice, and the optical properties are calculated using Mie theory (Whicker et al. 2022). Alternatively, the grains can be represented as spheres using Mie theory, as spheroids, hexagonal plates or koch snowflakes using calculations adapted from Mie theory (He et al. 2017, 2018), or hexagonal prisms using geometric optics calculations (van Diedenhoven 2014, Cook et al. 2020). The grain option works with both solvers as long as there is no fresnel layer included when using the AD solver. There is also an option to add liquid water coating around the grains using two-layer coated sphere Mie calculations (Niklas Bohn, 01/05/2020). This option had previously been included in the Matlab version of BioSNICAR and was deprecated for the Python translation in favour of increasing ice grain radii to simulate interstitial melt water; however, it was subsequently decided that giving the choice of method for incorporating liquid water was a better approach. Note that running the model with spherical grains and omitting any biological particles is equivalent to running the original SNICAR model of Flanner et al. (2007, 2009). The impact of colored dissolved organic matter (CDOM) on snow/ice optical properties can be incorporated by toggling on the cdom_layer parameter in the model. The current parametrisation is based on CDOM absorption measurements carried on samples from south-eastern greenlandic glaciers (Halbach et al. 2022).

### LAPs

The single scattering properties (ssps) of LAPs are calculated outside of the main driver of the model. For most of them, the ssps were calculated in published studies and are directly loaded from the main driver, so that the user only needs to select the type of LAP to include in the simulations from the driver (black carbon from Bohren and Huffman, 1983; brown carbon from Kirchstetter et al. 2004; dust from Balkanski et al 2007; volcanic ashes from Flanner et al 2014; Colorado dust from Skiles et al 2017; Greenlandic dust with low to high hematite content from Polashenski et al 2015, Greenlandic dust sampled in the south-western ablation area of the ice sheet from Cook et al. 2020; snow and glacier algae from Cook et al. 2017, Cook et al. 2020). 
However, BioSNICAR includes a bio-optical model that allows for the calculations of ssps of pigmented algae. Two options are available: spherical cells as per ice grains using Mie theory, typically used for snow algae, and circular based cylinders using geometric optics typically used for glacier algae that are long chains of cells approximated as cylinders after Lee and Pilon, 2013. The user needs to feed the bio-optical model with the refractive index and absorption cross section of the cell, along with algae dry and wet density and cell size. The absorption cross section can be calculated in the model from intracellular pigment concentrations (Cook et al. 2017; pigment absorption coefficients from Bidigare 1990, Halbach 2022, Clementson and Wojtasiewicz 2019)) or loaded as the absorbance of extracted pigments directly, which allows to use a pigment packaging correction (Chevrollier et al. 2022) or loaded as an empirical measurement of whole-cells. The ssps calculated are then stored into a file directly usable in the main driver of the model.

## Model Structure

<img src="./Assets/model_structure2.jpg" width=1500>


## Repository Structure

The following directory tree shows the correct structure for this model code. This is how the files are structured when this repository is cloned or downloaded. This can be used as a reference for understanding the software or as a guide if things get muddled during modification of the source code.

```
├── assets
│   ├── example-output.jpg
│   ├── model_schematic.odp
│   ├── model_structure2.jpg
│   ├── model_structure.jpg
│   ├── py_mat_comparison.png
│   └── SSA_derivation.pdf
|
├── bio_optical_model
│   ├── biooptical_driver.py
│   ├── biooptical_Funcs.py
│   ├── __pycache__
│   └── update_netCDFS.py
|
├── Data
│   ├── additional_data
│   │   ├── Albedo_master.csv
│   │   ├── ARF_master.csv
│   │   ├── Chlorophyll_a_m2mg.csv
│   │   ├── HCRF_master_16171819.csv
│   │   ├── InVivoPhenolData.csv
│   │   ├── PhenolicPigment_m2mg.csv
│   │   ├── phenol_mac_correction.csv
│   │   ├── phenol_MAC.csv
│   │   ├── phenol_mac_packaging_corrected.csv
│   │   ├── Photoprotective_carotenoids_m2mg.csv
│   │   ├── Photosynthetic_carotenoids_m2mg.csv
│   │   ├── pigmentMAC_200nm.csv
│   │   ├── pigmentMAC_250nm.csv
│   │   └── Spectra_Metadata.csv
│   ├── luts
│   │   └── LUT.npy
│   ├── OP_data
│   │   ├── 480band
│   │   │   ├── bubbly_ice_files
│   │   │   ├── fsds
│   │   │   ├── ice_hexagonal_columns
│   │   │   │   ├── ice_Pic16
│   │   │   │   ├── ice_Wrn08
│   │   │   │   └── ice_Wrn84
│   │   │   ├── ice_spherical_grains
│   │   │   │   ├── ice_Pic16
│   │   │   │   ├── ice_Wrn08
│   │   │   │   └── ice_Wrn84
│   │   │   ├── lap
│   │   │   ├── rain_polished_ice_spectrum.csv
│   │   │   └── r_sfc
│   │   │       ├── blue_ice_spectrum_s10290721.csv
│   │   │       └── rain_polished_ice_spectrum.csv
│   │   ├── ice_k.csv
│   │   ├── ice_n.csv
│   │   ├── ice_optical_constants.csv
│   │   ├── k_cdom_240_750.csv
│   │   ├── k_ice_480.csv
│   │   ├── Refractive_Index_Ice_Warren_1984.csv
│   │   ├── Refractive_Index_Liquid_Water_Segelstein_1981.csv
│   │   ├── water_RI.csv
│   │   └── wavelengths.csv
│   └── pigments
│       ├── alloxanthin.csv
│       ├── antheraxanthin.csv
│       ├── chl-a.csv
│       ├── chl-b.csv
│       ├── cis_astaxanthin_diester.csv
│       ├── cis_astaxanthin_monoester.csv
│       ├── lutein.csv
│       ├── neoxanthin.csv
│       ├── pckg_GA.csv
│       ├── pckg_SA.csv
│       ├── pheophytin.csv
│       ├── Photop_carotenoids.csv
│       ├── Photos_carotenoids.csv
│       ├── ppg.csv
│       ├── total_astaxanthin.csv
│       ├── trans_astaxanthin.csv
│       ├── trans_astaxanthin_ester.csv
│       ├── violaxanthin.csv
│       └── zeaxanthin.csv
├── experiments
│   ├── call_snicar.py
│   ├── config.yaml
│   ├── driver.py
│   ├── __init__.py
│   ├── __pycache__
│   ├── README.md
│   └── utils.py
|
├── LICENSE
├── README.md
├── requirements.txt
├── src
│   ├── adding_doubling_solver.py
│   ├── bubble_reff_calculator.py
│   ├── classes.py
│   ├── column_OPs.py
│   ├── display.py
│   ├── geometric_optics_ice.py
│   ├── __init__.py
│   ├── inputs.yaml
│   ├── mie_coated_water_spheres.py
│   ├── __pycache__
│   ├── setup_snicar.py
│   ├── snicar_driver.py
│   ├── toon_rt_solver.py
│   └── validate_inputs.py
└── tests
    ├── conftest.py
    ├── __init__.py
    ├── matlab_benchmark_script.m
    ├── __pycache__
    ├── README.md
    ├── snicarAD_v4.m
    ├── test_data
    │   ├── matlab_benchmark_data_clean.csv
    │   ├── matlab_benchmark_data.csv
    │   ├── py_benchmark_data_clean.csv
    │   ├── py_benchmark_data.csv
    │   └── py_mat_comparison.png
    └── test_snicar.py
```


# Testing

This repository contains a set of executable tests that anyone can run independently to verify the codebase against a Matlab benchmark version (published in Whicker et al 2021) along with a fuzzer that tests that the model runs and returns valid results across a wide parameter space. The fuzzer can be configured for specific ranges of values by adjusting the `pytest.mark.parameterize` statements in `test_snicar.py`, or just use our recommended defaults. The fuzzer can be toggled off by setting `fuzz = 0` in `conftest.py`. Tests are organised in `/tests` but are run from the top level directory so that the model code can be more conveniently imported as modules into the test session. Therefore, to run the tests, simply navigate to the top level directory and run:

`$ pytest tests`

This will open two datasets containing 5000 simulations replicated in the Python and Matlab implementations. Sucessfully passing tests are reported in the console as green dots, and pytest will return a summary of N tests passed and N tests failed. A figure showing N pairs of spectra is saved to the /tests folder for visual inspection. Any failures will be documented in the terminal so that they can be analysed and any bugs fixed. In the downloaded version of this repo, 100% of 1060 individual tests pass.

This demonstrates physically realistic predictions and equivalency between the two codebases to at least 1e-8 albedo units. The great majority of the simulations match to within 1e-12.

<img src="./Assets/py_mat_comparison.png" width=500>

More tests can and will be added over time (please feel free to contribute tests)!

The model configuration used to generate the data used to drive the automated tests can be found in the `matlab_benchmark_script.m` and `python_benchmark_script.py` files. The Python version calls functions in `py_benchmarking_funcs.py` The matlab version was run on a linux server at UMich, the Python version was run locally by the repository owner on Ubuntu 20.04. 

# Contributions

New issues and pull requests are welcomed. Pull requests trigger our Github Actions workflow to test for any breaking changes. PRs that pass these automated tests will be reviewed.

# Permissions

This code is provided under an MIT license with the caveat that it is in active development. Collaboration ideas and pull-requests generally welcomed. Please use the citations below to credit the builders of this repository and its predecessors.

# Citation

If you use this code in a publication, please cite:

Cook, J. et al. (2020): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere, doi:10.5194/tc-14-309-2020 

Flanner, M. et al. (2007): Present-day climate forcing and response from black carbon in snow, J. Geophys. Res., 112, D11202, https://doi.org/10.1029/2006JD008003

And if using the adding-doubling method please also cite Dang et al (2019) and Whicker et al (2022) as their code was translated to form the adding_doubling_solver.py script here. The aspherical grain correction equations comes from He et al. (2016).


# References

Balkanski, Y., Schulz, M., Claquin, T., & Guibert, S. (2007). Reevaluation of Mineral aerosol radiative forcings suggests a better agreement with satellite and AERONET data. Atmospheric Chemistry and Physics, 7(1), 81-95.

Bidigare, R. R., Ondrusek, M. E., Morrow, J. H., & Kiefer, D. A. (1990, September). In-vivo absorption properties of algal pigments. In Ocean Optics X (Vol. 1302, pp. 290-302). International Society for Optics and Photonics.

Bohren, C. F., & Huffman, D. R. (1983). Absorption and scattering of light by small particles. John Wiley & Sons.

Briegleb, B. P., and B. Light. "A Delta-Eddington multiple scattering parameterization for solar radiation in the sea ice component of the Community Climate System Model." NCAR technical note (2007).

Clementson, L. A., & Wojtasiewicz, B. (2019). Dataset on the absorption characteristics of extracted phytoplankton pigments. Data in brief, 24, 103875.

Cook JM, et al (2017) Quantifying bioalbedo: A new physically-based model and critique of empirical methods for characterizing biological influence on ice and snow albedo. The Cryosphere: 1–29. DOI: 10.5194/tc-2017-73, 2017b

Cook, J. M. et al. (2020): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-58, in review, 2019. 

Dang, C., Zender, C., Flanner M. 2019. Intercomparison and improvement of two-stream shortwave radiative transfer schemes in Earth system models for a unified treatment of cryospheric surfaces. The Cryosphere, 13, 2325–2343, https://doi.org/10.5194/tc-13-2325-2019 

Flanner, M. et al. (2007): Present-day climate forcing and response from black carbon in snow, J. Geophys. Res., 112, D11202, https://doi.org/10.1029/2006JD008003

Flanner, M et al. (2009) Springtime warming and reduced snow cover from
carbonaceous particles. Atmospheric Chemistry and Physics, 9: 2481-2497, 2009.

Flanner, M. G., Gardner, A. S., Eckhardt, S., Stohl, A., & Perket, J. (2014). Aerosol radiative forcing from the 2010 Eyjafjallajökull volcanic eruptions. Journal of Geophysical Research: Atmospheres, 119(15), 9481-9491.

Flanner, M. G. et al., SNICAR-ADv3: a community tool for modeling spectral snow albedo, Geosci. Model Dev., 14, 7673–7704, https://doi.org/10.5194/gmd-14-7673-2021, 2021.

He, C., Takano, Y., Liou, K. N., Yang, P., Li, Q., & Chen, F. (2017). Impact of snow grain shape and black carbon–snow internal mixing on snow optical properties: Parameterizations for climate models. Journal of Climate, 30(24), 10019-10036.

He, C., Liou, K. N., Takano, Y., Yang, P., Qi, L., & Chen, F. (2018). Impact of grain shape and multiple black carbon internal mixing on snow albedo: Parameterization and radiative effect analysis. Journal of Geophysical Research: Atmospheres, 123(2), 1253-1268.

Kirchstetter, T. W., Novakov, T., & Hobbs, P. V. (2004). Evidence that the spectral dependence of light absorption by aerosols is affected by organic carbon. Journal of Geophysical Research: Atmospheres, 109(D21).

Lee, E., & Pilon, L. (2013). Absorption and scattering by long and randomly oriented linear chains of spheres. JOSA A, 30(9), 1892-1900.

Picard, G., Libois, Q., & Arnaud, L. (2016). Refinement of the ice absorption spectrum in the visible using radiance profile measurements in Antarctic snow. The Cryosphere, 10(6), 2655-2672.

Polashenski et al. (2015): Neither dust nor black carbon causing apparent albedo decline in Greenland's dry snow zone: Implications for MODIS C5 surface reflectance, Geophys. Res. Lett., 42, 9319– 9327, doi:10.1002/2015GL065912, 2015.

Skiles, S. M., Painter, T., & Okin, G. S. (2017). A method to retrieve the spectral complex refractive index and single scattering optical properties of dust deposited in mountain snow. Journal of Glaciology, 63(237), 133-147.

Toon, O. B., McKay, C. P., Ackerman, T. P., and Santhanam, K. (1989), Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres, J. Geophys. Res., 94( D13), 16287– 16301, doi:10.1029/JD094iD13p16287. 

van Diedenhoven et al. (2014): A flexible paramaterization for shortwave opticalproperties of ice crystals. Journal of the Atmospheric Sciences, 71: 1763 – 1782, doi:10.1175/JAS-D-13-0205.1

Warren, S. G. (1984). Optical constants of ice from the ultraviolet to the microwave. Applied optics, 23(8), 1206-1225.

Warren, S. G., & Brandt, R. E. (2008). Optical constants of ice from the ultraviolet to the microwave: A revised compilation. Journal of Geophysical Research: Atmospheres, 113(D14).

Whicker et al., Halbach et al., Chevrollier et al. coming soon!
