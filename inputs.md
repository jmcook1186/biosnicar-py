# Input Configuration Guide for BioSNICAR

BioSNICAR is configured by updating the values in the config file `inputs.yaml`. This file is divided into subsections for variables that corrspond to specific parts of the model. These are:

1) `CTRL`: variables that control how the program executes ("metaconfig")
2) `RTM`: variables that configure the radiative transfer model
3) `ICE`: variables that adjust the ice physical properties
4) `PATHS`: paths to specific files
5) `IMPURITIES`: properties of impurities included in the simulation
6) `PLOT`: configuration for plotting output data

The user sets the desired values in `inputs.yaml` and then there is no need to interact with any other script in the repository apart from simpyl runnign the driver `python ./src/snicar_driver.py`. However, the value chpices are quite nuanced and require some explanation. In this document, the meaning of each field in `inputs.yaml` is defined along with any constraints and defaults.

## CTRL
### SMOOTH:
#### Definition
    Boolean that toggles smoothing of the predicted albedo. If `SMOOTH==True` a Savistky-Golay smoothign function is applied to the spectral albedo, with a window size and polynomial order defined by `WINDOW_SIZE` and `POLY_ORDER` respectively.

#### Values:
    `True` or `False`

#### Default:
    `False`

### WINDOW_SIZE:
#### Definition
    Interpolation window width used in Savisky-Golay filter if `SMOOTH==True`.

#### Values:
    Positive integer. Technically valid between 2-480 but practical between ~5-20.
#### Default:
    9


### POLY_ORDER:
#### Definition
   Order of polynomial interpolation function used in Savisky-Golay filter if `SMOOTH==True`.

#### Values:
    Positive integer. Technically any valid but practical between ~5-20.
#### Default
    3

## RTM
### DIRECT:
#### Definition
   Toggles between incoming irradiance being diffuse or direct. It's effect is to change a string stub that generates the filepath used to load in iorradiance from file. If direct irradiance is toggled the model reads in a file whose name contains "_clr_" and a solar zenith angle. 

#### Values:
    '0' for diffuse irradiance. `1` for direwct irradiance.
#### Default
    1

### NBR_WVL:
#### Definition
    Number of individual wavelengths considered by the model. Under normal circumstances the wavelength range is from 0.205 um to 4.995 um in steps of 10 nm, giving 480 individual wavelengths.

#### Values:
    positive integer
#### Default
    480

### APRX_TYP:
#### Definition
    Chooses two-stream approximation type. There are three choices: Eddington, Quadrature and Hemispheric mean. Each is a way to integrate over angle to give a single upwards and downwards flux. The hemispheric mean is derived by assuming the phase function to equal `1 + g` (where `g`= asymmetry parameter) in the forward scattering hemisphere and `1-g` in the backward scattering hemisphere. The hemispheric mean is useful for resolving instabilities in the infrared wavelengths. The equations for each of these approximations is provided in Toon et al. (1989) Table 1. 
    
    These choices are only available for the Toon solver - the adding-doubling solver uses hemispheric mean automatically.

#### Values:
    0: Eddington
    1: Quadrature
    2: Hemispheric mean

#### Default
    2

### DELTA:
#### Definition
    Toggles whether or not the delta function approximation is applied. The delta function truncates the forward scattering peak in the phase function and compensates by adjusting the voluem scattering coefficient, as described [here](https://ui.adsabs.harvard.edu/abs/1970JAtS...27..943P/abstract). 

    This choice is only available for the Toon solver - the adding-doubling applies the delta function automatically.

#### Values:
    True: Delta function
    False: No delta function
    
#### Default
    True


### SOLZEN:
#### Definition
    Solar zenoth angle. This refers to the angular position of the solar disc relative to vertical (nadir = 0, horizon = 90) expressed in degrees.

#### Values:
    positive integer degrees.
    For the Ton solver there are known instabilities for solzen outside the range 50-89. For the AD solver any positive integer between 1-89 is valid.
    
#### Default
    50

### INCOMING:
#### Definition
    Chooses spectal profile for incoming irradiance from one of seven presets. These are:
    
    0 = mid-latitude winter
    1 = mid-latitude summer
    2 = sub-Arctic winter
    3 = sub-Arctic summer
    4 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
    5 = High Mountain (summer, surface pressure of 556 hPa)
    6 = Top-of-atmosphere
Note that clear-sky spectral fluxes are loaded when direct_beam=1,
and cloudy-sky spectral fluxes are loaded when direct_beam=0

#### Values:
    Positive integer between 0 and 6
    
#### Default
    4






  DELTA: True
  SOLZEN: 50
  INCOMING: 0
  ILLUMINATION_FILE_STUBS: ["swnb_480bnd_mlw", "swnb_480bnd_mls", "swnb_480bnd_saw",
                          "swnb_480bnd_sas",  "swnb_480bnd_smm", "swnb_480bnd_hmn",
                          "swnb_480bnd_trp"]
  VIS_MAX_IDX: 39
  NIR_MAX_IDX: 480
