# BioSNICAR_GO_PY

Python implementation of BioSNICAR_GO, designed to predict the spectral and broadband albedo and energy absorption in snow and ice with variable grain sizes and densities either clean or with both inorganic and biological particles.


# Background

This is a Python implementation of the BioSNICAR_GO model offering the user a choice between determining ice optical properties using Mie scattering (good for fine snow grains that can be assumed spherical) or geometric optics (good for wet snow and ice). This also enables different grain shapes to be used - geometric optics for hexagonal plates and columns, Mie scattering for spheres, spheroids, hexagonal plates and Koch snowflakes (after He et al. 2016, 2017). This builds upon legacy FORTRAN and Matlab code developed by Flanner et al. (2007) and Cook et al. (2017, 2020). Note that running the model in Mie mode with spherical grains and omitting any biological particles is equivalent to running the original SNICAR model of Flanner et al. (2007, 2009).

The model is a two stream radiative transfer model that predicts the albedo of a snow or ice column with user-defined grain size and density and mass mixing ratios of a range of biological and non-biological particles. In GO mode, the optical properties of the ice grains themselves are calculated using a parameterisation of geometric optics calculations described by van Diedenhoven (2014) and the refractive indices are from Warren and Brandt (2008). In Mie mode the ice optical properties are calculated using Mie scatering codes using the same ice refractive indices. The rationale behind providing GO functionality is that geometric optics enables ice grains shaped as aritrarily large hexagonal plates and columns to be simulated, whereas Mie scattering applies to small spheres. The geometric optics approach is therefore better suited to wet snow or glacier ice, whereas Mie scattering is better for dry snow. Since 01/05/2020 there exists an option to model the effects of liquid water coatings around ice grains using two-layer coated sphere Mie calculations (added by Niklas Bohn). This option had previously been included in the Matlab version of BioSNICAR and was deprecated for the Python translation in favour of increasing ice grain radii to simulate interstitial melt water; however, it was subsequently decided that giving the choice of method for incorporating liquid water was a better approach. Gratitude to Niklas for adding this feature. 

BioSNICAR_GO also offers the option to incorporate algae as a light absorbing impurity, importantly either in the form of snow algae modelled using mie theory, or glacier algae modelled as large (10s - 100s microns), long chains of cells approximated as cylinders after Lee and Pilon, 2013) modelled using geometrical optics. 

The single scattering optical properties of the algae is determined by a separate Bio-optical model also provided here. Depending whether the user provides a MAC or pigment profile, the bio-optical model returns a cell MAC and imaginary refractive index that is parsed to the Mie or GO solver that in turn calculates the single scattering properties and appends them to a lookup library accessible to the radiative transfer model. The mass absorption coefficients for glacier algae provided in this repository have been determined empirically (Williamson et al. in review).

The mineral dusts included in this version include four "global average" dusts with typically Saharan optical properties as described by Flanner et al. 2009. There are three Greenland Ice Sheet mineral dusts that were generated from field measurements of local mineralogy on the south-western sector of the ice sheet near Kangerlussuaq. The optical properties for each mineral was obtained from the literature and mixed according to the measured relative abundance using the Maxwell-Garnett mixing model, and the optical properties predicted using Mie theory over a measured particle size distribution. There are also three additional Greenlandic dusts from Polashenski et al. (2015) who generated hypothetical dust samples with low, medium and high hematitie content, with the remainder being similar to dusts collected on Greenlandic snow. Black carbon is included in both sulphate-coated and uncoated forms, and volcanic ash is included as per Flanner et al. (2009).

This functionality, combined with the geometric optics option for the ice matrix makes BioSNICAR_GO applicable to bare glacier ice surfaces as well as wet and dry snowpacks - a significant step forwards from the original BiOSNICAR model published in Cook et al (2017). The model currently includes options for 2 glacier algae of different user-defined dimensions and one snow algae.

## Recent Additions

In 2020 there have been several additions to the BioSNICAR model that have not yet been documented in any publication paper. These include:

### 1) Refactor of bio-optical model
The bio-optical model was refactored into a much more useable format. Hard coded variables were moved into the function call, the two bio-optical components (the mixing model and mie/go code) ware now controlled from a single driver script, and file paths were all synchronised.

### 2) Update of glacier algae optical properties
The new glacier algal optical properties account for intracellular protein attachment and packaging of pigments into specific regions of the cell rather than assuming uniform distribution through the cell.
 
### 3) Addition of liquid water films
Thanks to Niklas Bohn to incorporating liquid water films of user-defined thickness to the model. This is controlled by providing a value for r_water when running the model in Mie mode and the optical properties are then calculated using a coated-spheres model.
   
### 4) Addition of specular reflection function
Specular reflection from the upper surface, scaled using a user-defined smoothness value is also now available. This approximates the the "specular-eddington" approach of Warren and Mullen (1984: doi 10.1029/JD093iD07p08403) and Dadic et al. (2013: doi 0.1002/jgrf.20098, 2013). If specular reflection is toggled ON, the model assumes the existence of a smooth layer at the upper surface where reflection occurs according to Fresnel's equations, giving a value for spectral specular reflectance (R) from the upper surface. This value is then used to attenuate the incoming irradiance that propagates to the two-stream multiple scattering model for n vertical layers below. The specular component is then added to the upwelling flux from the two-stream model when the material albedo is calculated.
    
(i.e. I* is total incoming irradiance, R is Fresnel reflectance from upper surface, I* x R is the specular component, I* x (1-R) is the irradiance propagated to the two-stream model, giving F_up, F-up_tot = (F_up + I* x R)).  

**** <b>CAUTION!!</b> Note that internal reflection (specular reflection of upwelling diffuse radiation from the underside of the upper ice/air boundary) is not yet accounted for - this is the next most pressing item on my TODO list but likely requires replacing the tridiagonal matrix solver in SNICAR with a different procedure ***


# Model Structure

<img src="./Assets/model_structure.jpg" width=1500>


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

snicar8d_mie.py: main snicar function using mie optical properties
SNICAR_driver.py: driver script for BioSNICAR_GO package
snicar_mie_tests.py: script for running unit tests against matlab version
BioSNICAR_py.yaml: exported environment
Files in /UnitTests/SnicarX/ are all benchmark spectral albedos for unit testing as predicted by the Matlab version of SNICAR.

## Data

In this repository the requisite datasets are divided into four classes, each contained within their own subdirectory in "/Data". These are:

1) Mie files: NetCDF files containing optical properties of ice crystals predicted by Mie theory. Also contains inorganic LAP optical properties.
   
2) GO files: NetCDF files containing optical properties of ice crystals predicted by geometrical optics. Also contains inorganic LAP optical properties.

3) Algal Optics: NetCDF files containing optical properties for algae of varying length and radius as predicted by geometrical optics (glacier algae) or Mie theory (snow algae)

4) Cell/pigment MACs: csv files containing the mass absorption coefficients for individual pigments or cells used in the Bio-optical model. These files are currently sved "loose" in the Data directory.


## Repository Structure
```
BioSNICAR_GO_PY
|
|----- snicar8d_mie.py
|----- SNICAR_driver.py
|----- snicar_mie_tests.py
|----- BioSNICAR_py.yaml
|
|-----BioOptical_Model
|           |
|           |---Algae_GO.py
|           |---BioOptical_Model.py
|
|-----IceOptical_Model
|           |
|           |---Geometric_Optics_Ice.py
|
|
|-----Data
|       |
|       |---Algal_Optical_Props
|       |      |
|       |      |--all algal single scattering props in .nc
|       |      |--mineral dust, ash and BC optical props in .nc
|       |
|       |---GO_files
|       |      |
|       |      |--all ice single scattering props in .nc
|       |      |--mineral dust, ash and BC optical props in .nc
|       |
|       |---Mie_files
|              |
|              |--all ice single scattering props in .nc
|              |--mineral dust, ash and BC optical props in .n
|
|
|-----UnitTests
|            |
|            |----Snicar_mie
|            |        |
|            |        |----apprx_albedo.csv
|            |        |----coszen_albedo.csv
|            |        |----delta_albedo.csv
|            |        |----direct_albedo.csv
|            |        |----dz_albedo.csv
|            |        |----R_sfc_albedo.csv
|            |        |----rds_albedo.csv
|            |        |----rho_albedo.csv
|            |
|            |
|            |----Snicar_GO
|                  |
|                  |----apprx_albedo.csv
|                  |----coszen_albedo.csv
|                  |----delta_albedo.csv
|                  |----direct_albedo.csv
|                  |----dz_albedo.csv
|                  |----R_sfc_albedo.csv
|                  |----rds_albedo.csv
|                  |----rho_albedo.csv
|
|-----Assets
|       |
|       |--model_structure.jpg
|

```

# How to use

Using the BioSNICAR_GO package is straightforwards. The user should not need to engage with any scripts apart from the driver (SNICAR_driver.py). In this script, the user defines values for all relevant variables. In-script annotations provide fine detail about the nature of each variable. Running this script will then report the broadband albedo to the console and plot the spectral albedo to a new window (if running interactively in Ipython or Jupyter) or save the figure if the user toggles savefigs == True. It is strongly recommended to avoid modifying the snicar8d_mie.py and snicar8d_GO.py scripts.

BioOptical_model.py and Algae_GO.py are used to generate new optical properties for algae to include in BioSNICAR_GO. There should be sufficient in-script annotations to guide this process, although it is not suggested to apply purely theoretical optical properties in real scientific use-cases without some empirical "ground-truthing" as there are few published datasets for the MAC of snow and glacier algae to verify them against. In BioSNICAR_GO the glacier algal optical properties were empirically measured by Williamson et al (2020: PNAS) using data from the Black and Bloom project.

# Known bugs and gotchas

1) Diffuse + Eddington:
   
The combination of Eddington 2-stream approximation and diffuse incident radiation causes the albedo to go negative at wavelengths > ~1.4 um. Recommend using quadature or hemispheric mean approximations when diffuse
incident radiation is required.

2) SZA limits
   
The two stream aproximation seems to fall apart at high zenith angles (>~0.57). This is common to all versions of SNICAR and is explained in Toon et al. (1989).

3) Snow algae
   
While glacier algae MACs have been determined empirically, there is only a hypothetical snow algae with potentially realistic pigemnt concentrations derived from the literature. A more accurate, empirically-derived set of single scattering optical properties for real snow algae is needed.


# Unit Testing

Unit testing has been carried out by running the newly translated model and comparing the predicted albedo with that predicted by the Matlab version for identical variable values. The Matlab verson of SNICAR is itself benchmarked against the original FORTRAN version that is implemented in CLM and thoroughly tested. The unit testing undertaken here is archived in the foldr wdir/Unit_Tests/. To run the tests, simply run the snicar_mie_tests.py script. The script holds all variable values constant at a set of defaults defined in-script, apart from a single test variable. The test variable is given a set of values to iterate through. For each test case, the identical case was run using the Matlab version and the output stored in an external csv file. For each test condition, the code simply tests whether the Python and Matlab implementation predicts the same value per unit wavelength. The test is passed if all the predicted values agree within 0.0001. There are separate scripts for testing the Mie and Geometrical Optics implementations.

This Python code was developed and run in Python 3.6.8 64-bit downloaded as part of Anaconda 4.7.11 and run on a Linux (Ubuntu 16.04 LTS) OS using VS Code 1.39.2 (and also repeated in PyCharm 2019.2.3). These are the only conditions under which the code has been tested. Under those conditions the unit tests return the following:

```

************ UNIT TESTS ***************
***************  MIE ******************
*****************************************

All tests run with default variable values apart from single 
test variable

Values benchmarked against Matlab version, tests pass if Matlab 
and Python versions are equal to within 
1e-4

************************

APPRX = 1 TEST PASSED WITH TOLERANCE 1e-4
APPRX = 2 TEST PASSED WITH TOLERANCE 1e-4
APPRX = 3 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

NB DIRECT TESTS RUN WITH APRX = 2 DUE TO KNOWN BUG WITH EDDINGTON + DIFFUSE

DIRECT = 0 TEST PASSED WITH TOLERANCE 1e-4
DIRECT = 1 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

DELTA = 0 TEST PASSED WITH TOLERANCE 1e-4
DELTA = 1 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

COSZEN = 0.3 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.35 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.4 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.45 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.5 TEST PASSED WITH TOLERANCE 1e-4

 **************************************************

DZ = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #4 TEST PASSED WITH TOLERANCE 1e-4

 **************************************************

R_sfc = 0.2 TEST PASSED WITH TOLERANCE 1e-4
R_sfc = 0.4 TEST PASSED WITH TOLERANCE 1e-4
R_sfc = 0.6 TEST PASSED WITH TOLERANCE 1e-4

 **************************************************

RDS = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
RDS = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
RDS = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
RDS = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4
RDS = ARRAY #4 TEST FAILED WITH TOLERANCE 1e-4

 **************************************************

RHO = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #4 TEST PASSED WITH TOLERANCE 1e-4

**************************************************
```

```
************ UNIT TESTS ***************
********* Geometric Optics ************
*****************************************

All tests run with default variable values apart from single 
test variable

Values benchmarked against Matlab version, tests pass if Matlab 
and Python versions are equal to within 
1e-4

************************

NB DIRECT TESTS RUN WITH APRX = 2 DUE TO KNOWN BUG WITH EDDINGTON + DIFFUSE

DIRECT = 0 TEST PASSED WITH TOLERANCE 1e-4
DIRECT = 1 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

APPRX = 1 TEST PASSED WITH TOLERANCE 1e-4
APPRX = 2 TEST PASSED WITH TOLERANCE 1e-4
APPRX = 3 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

DELTA = 0 TEST PASSED WITH TOLERANCE 1e-4
DELTA = 1 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

COSZEN = 0.3 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.35 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.4 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.45 TEST PASSED WITH TOLERANCE 1e-4
COSZEN = 0.5 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

DZ = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4
DZ = ARRAY #4 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

R_sfc = 0.2 TEST PASSED WITH TOLERANCE 1e-4
R_sfc = 0.4 TEST PASSED WITH TOLERANCE 1e-4
R_sfc = 0.6 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

SIDE LENGTH = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
SIDE LENGTH = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
SIDE LENGTH = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
SIDE LENGTH = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

DEPTH = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
DEPTH = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
DEPTH = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

RHO = ARRAY #0 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #1 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #2 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #3 TEST PASSED WITH TOLERANCE 1e-4
RHO = ARRAY #4 TEST PASSED WITH TOLERANCE 1e-4

**************************************************

```

It would be sensible to check that the same results are obtained on a new local machine before using BioSNICAR_GO in earnest.

# Permissions

This code is in active development. Collaboration ideas and pull-requests generally welcomed. The author assumes no responsibility for downstream usage and, while collaboration is welcome, there is no obligation to provide training or technical support for users. This code is provided without warranty or guarantee, nor is there the implied warranty of merchantability or fitness for any particular purpose.

# Citation

Please cite:

Cook, J. et al. (2019): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere, https://doi.org/10.5194/tc-2019-58, accepted for publication Dec 2019. 

and the doi for the relevant release of this repository (v0.1 doi: 10.5281/zenodo.3564517).

# References

Cook JM, et al (2017) Quantifying bioalbedo: A new physically-based model and critique of empirical methods for characterizing biological influence on ice and snow albedo. The Cryosphere: 1–29. DOI: 10.5194/tc-2017-73, 2017b

Cook, J. M. et al. (2019): Glacier algae accelerate melt rates on the western Greenland Ice Sheet, The Cryosphere Discuss., https://doi.org/10.5194/tc-2019-58, in review, 2019. 

Flanner, M. et al. (2007): Present-day climate forcing and response from black carbonin snow, J. Geophys. Res., 112, D11202, https://doi.org/10.1029/2006JD008003

Flanner, M et al. (2009) Springtime warming and reduced snow cover from
carbonaceous particles. Atmospheric Chemistry and Physics, 9: 2481-2497, 2009.

Polashenski et al. (2015): Neither dust nor black carbon causing apparent albedo decline in Greenland's dry snow zone: Implications for MODIS C5 surface reflectance, Geophys. Res. Lett., 42, 9319– 9327, doi:10.1002/2015GL065912, 2015.

van Diedenhoven et al. (2014): A flexible paramaterization for shortwave opticalproperties of ice crystals. Journal of the Atmospheric Sciences, 71: 1763 – 1782, doi:10.1175/JAS-D-13-0205.1
