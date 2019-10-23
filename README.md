# BioSNICAR_GO_PY

Translation of the BioSNICAR_GO model into Python (including translation of the original SNICAR model)

# Current Development status

21 Oct 2019: Mie scattering version of snicar functional. No unit testing done yet, pending access to previous SNICAR versions for benchmarking. Geometrical optics versions are not yet translated into Python and therefore setting GeometricOptics = 1 in the driver script raises an exception and will not run.

# In this repo

## Scripts
snicar8d_mie.py: main snicar function using mie optical properties
SNICAR_driver.py: driver script for BioSNICAR_GO package
snicar_mie_tests.py: script for running unit tests against matlab version

## Folders
Unit_Tests: contains csv files for albedo predicted by Matlab version as target values for unit testing

## Repository Structure

BioSNICAR_GO_PY
|
|----- snicar8d_mie.py
|----- SNICAR_driver.py
|----- snicar_mie_tests.py
|
¬¬¬UnitTests
      |
      |----Snicar_mie
              |
              |----apprx_albedo.csv
              |----coszen_albedo.csv
              |----delta_albedo.csv
              |----direct_albedo.csv
              |----dz_albedo.csv
              |----R_sfc_albedo.csv
              |----rds_albedo.csv
              |----rho_albedo.csv


# Background

# How to use

# Known bugs and gotchas

1) Diffuse + Eddington
The combination of Eddington 2-stream approximation and diffuse incident radiation causes the albedo to go negative at wavelengths > ~1.4 um. Recommend using quadature or hemispheric mean approximations when diffuse
incident radiation is required.

2) SZA limits
The two stream aproximation seems to fall apart at high zenith angles (>~0.57). This is common to all versions of SNICAR and is explained in Toon et al. (1989).

3) Snow algae
While glacier algae MACs have been determined empirically, we have included only a hypothetical snow algae with potentially realistic pigemnt concentrations derived from the literature. A more accurate, empirically-derived set of single scattering optical properties for real snow algae is needed.


# Unit Testing

Unit testing has been carried out by running the newly translated model and comparing the predicted albedo with that predicted by the Matlab version for identical variable values. The Matlab verson of SNICAR is itself benchmarked against the original FORTRAN version that is implemented in CLM and thoroughly tested. The unit testing undertaken here is archived in the foldr wdir/Unit_Tests/. To run the tests, simply run the snicar_mie_tests.py script. The script holds all variable values constant at a set of defaults defined in-script, apart from a single test variable. The test variable is given a set of values to iterate through. For each test case, the identical case was run using the Matlab version and the output stored in an external csv file. For each test condition, the code simply tests whether the Python and Matlab implementation predicts the same value per unit wavelength. The test is passed if all the predicted values agree within 0.0001. 

This Python code was developed and run in Python 3.6.8 64-bit downloaded as part of Anaconda 4.7.11 and run on a Linux (Ubuntu 16.04 LTS) OS using VS Code 1.39.2 (and also repeated in PyCharm 2019.2.3). Under those conditions the unit test returns to following:

```

************ UNIT TESTS ***************

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

It would be sensible to check that the same results are obtained on a new local machine.

# Permissions
This code is in active development and no permissions are granted at this time.


# Citation
A release tag with a doi will be minted when this code has been fully tested, until then I do not recommend using this code.