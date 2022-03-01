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
#### Definition:
Boolean that toggles smoothing of the predicted albedo. If `SMOOTH==True` a Savistky-Golay smoothign function is applied to the spectral albedo, with a window size and polynomial order defined by `WINDOW_SIZE` and `POLY_ORDER` respectively.

#### Values:
`True` or `False`

#### Default:
`False`

### WINDOW_SIZE:
#### Definition:
Interpolation window width used in Savisky-Golay filter if `SMOOTH==True`.

#### Values:
Positive integer. Technically valid between 2-480 but practical between ~5-20.
#### Default:
9


### POLY_ORDER:
#### Definition:
Order of polynomial interpolation function used in Savisky-Golay filter if `SMOOTH==True`.

#### Values:
Positive integer. Technically any valid but practical between ~5-20.
#### Default:
3

## RTM
### DIRECT:
#### Definition:
Toggles between incoming irradiance being diffuse or direct. It's effect is to change a string stub that generates the filepath used to load in iorradiance from file. If direct irradiance is toggled the model reads in a file whose name contains "_clr_" and a solar zenith angle. 

#### Values:
'0' for diffuse irradiance. `1` for direwct irradiance.
#### Default:
1

### NBR_WVL:
#### Definition
Number of individual wavelengths considered by the model. Under normal circumstances the wavelength range is from 0.205 um to 4.995 um in steps of 10 nm, giving 480 individual wavelengths.

#### Values:
positive integer
#### Default:
480

### APRX_TYP:
#### Definition:
Chooses two-stream approximation type. There are three choices: Eddington, Quadrature and Hemispheric mean. Each is a way to integrate over angle to give a single upwards and downwards flux. The hemispheric mean is derived by assuming the phase function to equal `1 + g` (where `g`= asymmetry parameter) in the forward scattering hemisphere and `1-g` in the backward scattering hemisphere. The hemispheric mean is useful for resolving instabilities in the infrared wavelengths. The equations for each of these approximations is provided in Toon et al. (1989) Table 1. 
    
These choices are only available for the Toon solver - the adding-doubling solver uses hemispheric mean automatically.

#### Values:
0: Eddington
1: Quadrature
2: Hemispheric mean

#### Default:
2

### DELTA:
#### Definition:
Toggles whether or not the delta function approximation is applied. The delta function truncates the forward scattering peak in the phase function and compensates by adjusting the voluem scattering coefficient, as described [here](https://ui.adsabs.harvard.edu/abs/1970JAtS...27..943P/abstract). 

This choice is only available for the Toon solver - the adding-doubling applies the delta function automatically.

#### Values:
True: Delta function
False: No delta function
    
#### Default:
True


### SOLZEN:
#### Definition:
Solar zenith angle. This refers to the angular position of the solar disc relative to vertical (nadir = 0, horizon = 90) expressed in degrees.

#### Values:
Positive integer degrees.
For the Ton solver there are known instabilities for solzen outside the range 50-89. For the AD solver any positive integer between 1-89 is valid.
    
#### Default:
50

### INCOMING:
#### Definition:
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
    
#### Default:
4


### ILLUMINATION_FILE_STUBS:
#### Definition:
These file stubs are used to select the correct file to open for the incoming irradiance depending on the value of `incoming` and `solzen`.

#### Values:
The file stubs are strings and they are constant - do not change.

### Default:
["swnb_480bnd_mlw", "swnb_480bnd_mls", "swnb_480bnd_saw", "swnb_480bnd_sas",  "swnb_480bnd_smm", "swnb_480bnd_hmn", "swnb_480bnd_trp"]


### VIS_MAX_IDX:
#### Definition:
The index in `wavelengths` of the upper wavelength in the visible range.

#### Values:
Technically any positive integer between 1 and 480 is valid, but practically values close to 39 are appropriate as this corresponds to `wavelength` == 0.7um.

### Default:
39 

### NIR_MAX_IDX:
#### Definition:
The index in `wavelengths` of the upper wavelength in the near-infrared range.

#### Values:
Technically any positive integer between 1 and 480 is valid, but practically values close to 480 are appropriate as this corresponds to `wavelength` == 4.995 um.

### Default:
480 


## ICE
### DZ:
#### Definition:
Thickness of each layer in the model in units of meters.

#### Values:
array with length == nbr_lyr. Each element is the thickness of a layer in meters. `dz[0]` is the top layer.

### Default:
[0.1, 0.1] 

### LAYER_TYPE
#### Definition:
There are two layer tupes available in the model: granular or solid ice. Granular ice corresponds to a bulk medium of air with discrete ice grains. Solid ice refers to ice being the bulk medium with air inclusions. Solid ice layers have Fresnel reflection at their upper surfaces.

The toon solver can only accept granular layers. The adding-doubling solver accepts either type including a mix of both in any order.

#### Values:
0 for granular ice, 1 for solid ice

### Default:
[0, 0] 


### RHO
#### Definition:
Ice density in each layer measured in kg/m3. 

#### Values:
Positive integers up to 917 (density of pure ice).

### Default:
[400, 400] 

### RF
#### Definition:
Choice of refractive index for ice. The refractive index of pure ice has been measured several times over the past 5 decades. Here we toggle between several of those datasets. The options are:

0: Warren et al. (1984)
1: Warren et al. (2008)
2: Picard et al. (2016)

#### Values:
Positive integer between 0 and 2.

### Default:
2 

### CDOM
#### Definition:
Chromophoric dissolved organic matter (CDOM) is coloured organic matter that can be found in glacier meltwater. This variable toggles its presence on or off in each layer. This is an experimental feature that only allows a single concentration and it can only be added to solid ice layers (layer_type=1). The CDOM absorption data come from an upcoming paper by Halbach et al. Toggling CDOM on means adjustign the imaginary refractive index of the ice for the presence of CDOM, assuming it to well mixed through the ice layer.

#### Values:
0 for no CDOM, 1 for CDOM

### Default:
[0, 0] 


### SHP
#### Definition:
Shape of ice grains. Adjustments for grain shape can onyl be made when layer_type ==0. The available shapes are:

0 = sphere
1 = spheroid
2 = hexagonal plate
3 = koch snowflake
4 = hexagonal prisms


#### Values:
Positive integer between 0 and 4

### Default:
[0, 0] 


  SHP: [0, 0]  # grain shape (He et al. 2016, 2017)
  RDS: [12500, 12500]  # effective grain or bubble radius
  WATER: [0, 0]  # radius of optional liquid water coating
  HEX_SIDE: [ 10000, 10000]
  HEX_LENGTH: [10000, 10000]
  SHP_FCTR: [ 0, 0]
  AR: [0, 0]