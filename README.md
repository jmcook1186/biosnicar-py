# BioSNICAR_GO_PY 


<b>NOTICE (Jan 2022): We have done a lot of work refactoring and updating this code recently to make it cleaner, faster and more user-friendly. Some of the docs are currently lagging behind - we are working on this - please bear with us! </b>


This package is a Python implementation of the BioSNICAR_GO model to represent the albedo of snow and ice with variable grain sizes and densities, either clean or with light absorbing particulates  (LAPs). It builds upon legacy FORTRAN and Matlab code developed by Flanner et al. (2007), Cook et al. (2017, 2020) and Wicker et al. (2021, in review). The model is a two stream radiative transfer model that predicts the albedo of a snow or ice multilayer column with user-defined parameters for the snow/ice background (single scattering properties, density, grain size) as well as the LAPs (single scattering properties, concentration in the ice). Two solvers are available in this new version, the original SNICAR matrix solver typically representing ice and snow with grains (Toon et al. 1989) and the Adding-Doubling (AD) solver (Wicker et al. 2021, in review), representing the ice as a medium including air bubbles and allowing the incorporation of Fresnel layers. It also includes a BioOptical model specifically designed to represent biological albedo reducers on snow and ice surfaces. This functionality, combined with the AD solver for the ice matrix makes BioSNICAR_GO applicable to bare glacier ice surfaces as well as wet and dry snowpacks - a significant step forwards from the original BiOSNICAR model published in Cook et al (2017). 

# Background
## Snow and ice 

The snow and ice single scattering properties can be generated from the ice refractive index (Warren 1984, Warren and Brandt 2008, Picard 2016) using Mie scattering (good for fine snow grains that can be assumed spherical, or for representing ice bubbles in an ice matrix) or geometric optics (adapted from van Diedenhoven (2014), good for wet snow and ice). This also enables different grain shapes to be used - geometric optics for hexagonal plates and columns, Mie scattering for spheres, spheroids, hexagonal plates and Koch snowflakes (after He et al. 2016, 2017). The rationale behind providing GO functionality is that geometric optics enables ice grains shaped as large hexagonal plates and columns to be simulated, whereas Mie scattering applies to small spheres. Note that running the model in Mie mode with spherical grains and omitting any biological particles is equivalent to running the original SNICAR model of Flanner et al. (2007, 2009). Since 01/05/2020 there exists an option to model the effects of liquid water coatings around ice grains using two-layer coated sphere Mie calculations (added by Niklas Bohn). This option had previously been included in the Matlab version of BioSNICAR and was deprecated for the Python translation in favour of increasing ice grain radii to simulate interstitial melt water; however, it was subsequently decided that giving the choice of method for incorporating liquid water was a better approach. Gratitude to Niklas for adding this feature.  For each layer, the density has to be indicated together with the shape, coating and size of the grains (ice or air bubbles).

## LAPs

The LAPs single scattering properties can be also be generated using Mie theory (good for mineral dust or snow algae) or geometric optics (good for glacier algae that are long chains of cells approximated as cylinders after Lee and Pilon, 2013). In particular, the BioSNICAR_GO package includes a BioOptical model designed to generate the optical properties of snow and glacier algae cells and store them in a file directly usable in the main driver. The user needs to prescribe algae dry density, real part of refractive index and cell size for the calculation of the scattering properties, and for the calculation of extinction properties the cells are assumed homogeneous so that they are directly calculated from absorption cross sections (ACS), prescribed by an empirical measurement or a pigment profile (cellular pigment concentrations) from which the ACS is reconstructed following Cook et al. 2017, Pottier et al. 2005. The ACS of both glacier and snow algae are to date only theoretically reconstructed from pigment profiles.
The mineral dusts included in this version include four "global average" dusts with typically Saharan optical properties as described by Flanner et al. 2009. There are three Greenland Ice Sheet mineral dusts that were generated from field measurements of local mineralogy on the south-western sector of the ice sheet near Kangerlussuaq. The optical properties for each mineral was obtained from the literature and mixed according to the measured relative abundance using the Maxwell-Garnett mixing model, and the optical properties predicted using Mie theory over a measured particle size distribution. There are also three additional Greenlandic dusts from Polashenski et al. (2015) who generated hypothetical dust samples with low, medium and high hematitie content, with the remainder being similar to dusts collected on Greenlandic snow. Black carbon is included in both sulphate-coated and uncoated forms, and volcanic ash is included as per Flanner et al. (2009).

## Recent Additions

There have been several additions to the BioSNICAR model that have not yet been documented in any publication paper. These include:

### 1) Refactor of bio-optical model
The bio-optical model was refactored into a much more useable format. Hard coded variables were moved into the function call, the two bio-optical components (the mixing model and mie/go code) ware now controlled from a single driver script, and file paths were all synchronised.

### 2) Update of glacier algae optical properties
The new glacier algal optical properties account for intracellular protein attachment and packaging of pigments into specific regions of the cell rather than assuming uniform distribution through the cell using a linear correction.
 
### 3) Addition of liquid water films
Thanks to Niklas Bohn to incorporating liquid water films of user-defined thickness to the model. This is controlled by providing a value for r_water when running the model in Mie mode and the optical properties are then calculated using a coated-spheres model.

### 4) October 2020: Addition of aspherical grains
Adjustments to the single scattering optical properties of the ice grains resulting from variations in grain shape have been added to the model when run in Mie mode - running in GO mode assumes hexagonal columnar grains. These adjustments are from He et al. (2016) and enable the modelling of spheroids, hexagonal plates and Koch snowflakes.

### 5) November 2020: Addition of adding-doubling solver leading to solid ice layers and fresnel reflection
A second option for the treatment of multiple layers was added in October 2020. The user can now opt to include solid ice layers with fresnel reflections at their surfaces and solve the radiative transfer using the adding-doubling method as an alternative to the Toon et al (1989) matrix inversion method that was previously used that cannot include fresnel reflection and refraction. The A-D solver is pretty much a direct translation of the Matlab code written by Chloe Whicker (UMich) who in turn built on work by Dang et al. (2019). The solid ice layers can be positioned anywhere in the vertical profile, meaning unweathered ice beneath an increasingly porous weathered crust can now be realistically simulated.

### 6) November 2020: removed separate scripts for GO and Mie modes and synthesised into single SNICAR_feeder script, extended model to 200nm and enabled GO ice crytals in AD mode
Model now runs across 480 wavelengths between 200 - 5000 nm in all configurations, including GO.
The model now works using a single call to a snicar feeder script where the various radiative transfer properties (tau, SSA, g etc) are calculated and then passed to one of the two solvers.

### 7) January 2021: added testing scripts for comparing new code against matlab benchmarks. Extend spectral range of all model modes to 200nm
Reorganised directory structure to gather all testing scripts and data into the folder "Unit_Tests". Using Whicker/Flanner's Matlab code as benchmark for testing Python code.

### 8) September 2021: Bio-optical model renewed, new units
The bio-optical model was updated into a much simpler and user-friendly format that can be run from a short driver script. In addition, the user now has the option to generate or input the algae ACS in m2/cell or m2/um3 from the pigment profiles. This allows to reduce the uncertainty on the final ACS, as cell numbers and biovolumes are empirically determined, in contrary to single cell mass.

### 9) October 2021: CDOM layers and algae concentrations in cells/mL
The main driver now includes the possibility to indicate algae concentrations in cells/mL through the "GA/SA_units" variable, which is the standard unit for algae quantification from environmental samples. In this mode, the ACS needs to be indicated in m2/cell. Experimental CDOM layers have also been added based on coefficients from Halbach et al. 2021 (in prep), so that the impact of CDOM can be represented in each layer of the snow/ice column. 


# Model Structure

<img src="./Assets/model_structure2.jpg" width=1500>


# Environment/Dependencies

It is recommended to run this code in a fresh environment as follows:

```
conda create -n BioSNICAR_py python 3.6.8, numpy, scipy, matplotlib, pandas
conda install xarray dask netCDF4 bottleneck
pip install miepython

```

or alternatively linux users can create from the yaml file provided in this repository (which also includes packages such as jupyter for running interactively in vscode).

```
conda env create -f BioSNICAR_py.yaml

```
If creating a new environment in this way, please remember to update the prefix in the final line of the file to something appropriate to the new local machine. The .yaml file has been tested on Ubuntu 16.04 and 20.04.

# In this repo

## Files

SNICAR_driver.py: driver script for BioSNICAR_GO package
SNICAR_feeder.py: script called from SNICAR_driver, itself calling the solvers
Toon_RT_solver.py : Toon et al (1989) tridiagonal matrix solver
adding_doubling_solver.py: adding-doubling method radiative transfer solver
snicar_mie_tests.py: script for running unit tests against matlab version
BioSNICAR_py.yaml: exported environment

## Data

In this repository the requisite datasets are divided into four classes, each contained within their own subdirectory in "/Data". These are:

1) Mie files: NetCDF files containing optical properties of ice crystals predicted by Mie theory. Also contains all LAP optical properties.
   
2) GO files: NetCDF files containing optical properties of ice crystals predicted by geometrical optics. 

3) Algal_Optical_Props: NetCDF files containing optical properties for algae of varying length and radius as predicted by geometrical optics (glacier algae) or Mie theory (snow algae)

4) pigments: csv files containing the mass absorption coefficients for individual pigments or cells used in the Bio-optical model. 

5) bubbly_ice_files: optical properties for slabs of ice wth bubbles. numeric value in filename is the
   effective radius of the bulk material (higher number  = fewer bubbles)

## Repository Structure
```
BioSNICAR_GO_PY
|
|----- snicar_driver.py
|----- snicar_feeder.py
|----- Toon_RT_solver.py
|----- adding_doubling_solver.py
|----- snicar_mie_tests.py: script for running unit tests against matlab version
|----- BioSNICAR_py.yaml: exported environment
|
|-----BioOptical_Model
|           |
|           |---BioOptical_driver.py 
|           |---BioOptical_Funcs.py
|
|-----IceOptical_Model
|           |
|           |---Geometric_Optics_Ice.py
|           |---mie_coated_water_spheres.py
|
|
|-----Data
|       |
|       |- pigment MACs as csv files, ice and water refractive indices as .nc file
|       |
|       |---Algal_Optical_Props
|       |      |
|       |      |--algal single scattering props in .nc
|       |
|       |---GO_files
|       |      |
|       |      |--all ice single scattering props in .nc
|       |
|       |---Mie_files
|       |      |
|       |      |--all ice single scattering props in .nc
|       |      |--mineral dust, ash, BC, algae optical props in .nc
|       |
|       |---bubbly_ice_files
|              |
|              |--optical properties for ice slabs
|
|----------tests
|            |
|            IceOpticalModel (to match main repo)
|            snicar feeder.py (to match main repo)
|            adding_doubling_solver.py (to match main repo)
|            matlab_benchmark_script.m      
|            py_benchmark_script.py  
|            matlab_benchmark_data.csv
|            py_benchmark_data.csv
|            test_snicar.py
|            conftest.py
|            benchmarking_funcs.py
|            constants.txt
|            variables.txt
|
|-----Assets
|       |
|       |--model_structure.jpg
|       |--other images
|

```

# How to use

In top level directory:

`python ./src/snicar_driver.py`


Using the BioSNICAR_GO package is straightforwards. The user should not need to engage with any scripts apart from the driver (snicar_driver.py). In this script, the user defines values for all relevant variables. In-script annotations provide fine detail about the nature of each variable. Running this script will then report the broadband albedo to the console and plot the spectral albedo to a new window (if running interactively in Ipython or Jupyter) or save the figure if the user toggles savefigs == True. It is strongly recommended to avoid modifying the SNICAR_feeder script.

BioOptical_driver.py is used to generate new optical properties for algae to include in BioSNICAR_GO. There should be sufficient in-script annotations to guide this process, although it is not suggested to apply purely theoretical optical properties in real scientific use-cases without some empirical "ground-truthing" as there are few published datasets for the MAC of snow and glacier algae to verify them against. In BioSNICAR_GO the glacier algal optical properties were generated using empirically measured pigment profiles by Williamson et al (2020: PNAS) using data from the Black and Bloom project.


# Testing

This repository contains a set of executable tests that anyone can run independently to verify the codebase against a Matlab benchmark version (published in Whicker et al 2021) along with a fuzzer that tests that the model runs and returns valid results across a wide parameter space. The fuzzer can be configured for speicific ranges of values by adjusting the `pytest.mark.parameterize` statements in `test_snicar.py`, or just use our recommended defaults. The fuzzer can be toggled off by setting `fuzz = 0` in `conftest.py`. Tests are organised in `/tests` but are run from the top level directory so that the model code can be more conveniently imported as modules into the test session. Therefore, to run the tests, simply navigate to the top level directory and run:

`$ pytest tests`

This will open two datasets containing 5000 simulations replicated in the Python and Matlab implementations. Sucessfully passing tests are reported in the console as green dots, and pytest will return a summary of N tests passed and N tests failed. A figure showing N pairs of spectra is saved to the /tests folder for visual inspection. Any failures will be documented in the terminal so that they can be analysed and any bugs fixed. In the downloaded version of this repo, 100% of 1060 individual tests pass.

This demonstrates physically realistic predictions and equivalency between the two codebases to at least 1e-8 albedo units. The great majority of the simulations match to within 1e-12.

<img src="./Assets/py_mat_comparison.png" width=500>

More tests can and will be added over time (please feel free to contribute tests)!

The model configuration used to generate the data used to drive the automated tests can be found in the `matlab_benchmark_script.m` and `python_benchmark_script.py` files. The Python version calls functions in `py_benchmarking_funcs.py` The matlab version was run on a linux server at UMich, the Python version was run locally by the repository owner on Ubuntu 20.04. 


# Permissions

This code is in active development. Collaboration ideas and pull-requests generally welcomed. The author assumes no responsibility for downstream usage and, while collaboration is welcome, there is no obligation to provide training or technical support for users. This code is provided without warranty or guarantee, nor is there the implied warranty of merchantability or fitness for any particular purpose.

# Citation

If you use this code in a publication, please cite:

Cook, J. et al. (2020): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere, doi:10.5194/tc-14-309-2020 

Flanner, M. et al. (2007): Present-day climate forcing and response from black carbon in snow, J. Geophys. Res., 112, D11202, https://doi.org/10.1029/2006JD008003

And if using the adding-doubling method please also cite Dang et al (2019) and Whicker et al (2021) as their code was translated to form the adding_doubling_solver.py script here. The aspherical grain correction equations come from He et al. (2016).


# References

Cook JM, et al (2017) Quantifying bioalbedo: A new physically-based model and critique of empirical methods for characterizing biological influence on ice and snow albedo. The Cryosphere: 1–29. DOI: 10.5194/tc-2017-73, 2017b

Cook, J. M. et al. (2019): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-58, in review, 2019. 

Dang, C., Zender, C., Flanner M. 2019. Intercomparison and improvement of two-stream shortwave radiative transfer schemes in Earth system models for a unified treatment of cryospheric surfaces. The Cryosphere, 13, 2325–2343, https://doi.org/10.5194/tc-13-2325-2019 

Flanner, M. et al. (2007): Present-day climate forcing and response from black carbon in snow, J. Geophys. Res., 112, D11202, https://doi.org/10.1029/2006JD008003

Flanner, M et al. (2009) Springtime warming and reduced snow cover from
carbonaceous particles. Atmospheric Chemistry and Physics, 9: 2481-2497, 2009.

He, C., Liou, K.‐N., Takano, Y., Yang, P., Qi, L., & Chen, F. (2018). Impact of grain shape and multiple black carbon internal mixing on snow albedo: Parameterization and radiative effect analysis. Journal of Geophysical Research: Atmospheres, 123, 1253– 1268. https://doi.org/10.1002/2017JD027752 

Polashenski et al. (2015): Neither dust nor black carbon causing apparent albedo decline in Greenland's dry snow zone: Implications for MODIS C5 surface reflectance, Geophys. Res. Lett., 42, 9319– 9327, doi:10.1002/2015GL065912, 2015.

Toon, O. B., McKay, C. P., Ackerman, T. P., and Santhanam, K. (1989), Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres, J. Geophys. Res., 94( D13), 16287– 16301, doi:10.1029/JD094iD13p16287. 

van Diedenhoven et al. (2014): A flexible paramaterization for shortwave opticalproperties of ice crystals. Journal of the Atmospheric Sciences, 71: 1763 – 1782, doi:10.1175/JAS-D-13-0205.1

Whicker et al. COMING SOON!