#!/usr/bin/python

"""Classes used in BioSNICAR.

This module contains class definitions for all the classes used in BioSNICAR. This includes:

Impurity
Ice
Illumination
RTConfig
ModelConfig
PlotConfig
Outputs
BioOpticalConfig

These classes are used as convenient, mutable containers for the necessary data required to run
BioSNICAR. They are automatically instantiated by calling setup_snicar() using values provided
in inputs.yaml. Class functions are available for recalculating derived attributes when the
user changes attributes of Ice or Illumination classes. 

"""


import sys

sys.path.append("./src")
import math
import os

import numpy as np
import xarray as xr
import yaml


class Impurity:
    """Light absorbing impurity.

    Instances of Impurity are one discrete type of light absorbing impurity with
    a distinct set of optical properties.

    Attributes:
        name: name of impurity
        cfactor: concentration factor used to convert field measurements to model config (default=1)
        unit: the unit the concentration should be represented in (0 = ppb, 1 = cells/mL)
        conc: concentration of the impurity in each layer (in units of self.unit)
        file: name of netCDF file containing optical properties and size distribution
        impurity_properties: instance of opened file self.file
        mac: mass absorption coefficient (m2/kg or m2/cell)
        ssa: single scattering albedo
        g: asymmetry parameter

    """

    def __init__(self, file, coated, cfactor, unit, name, conc):

        self.name = name
        self.cfactor = cfactor
        self.unit = unit
        self.conc = conc
        self.file = file

        self.impurity_properties = xr.open_dataset(
            str(os.getcwd() + "/" + "/Data/OP_data/480band/lap/" + file)
        )

        if coated:
            mac_stub = "ext_cff_mss_ncl"
        else:
            mac_stub = "ext_cff_mss"

        self.mac = self.impurity_properties[mac_stub].values
        self.ssa = self.impurity_properties["ss_alb"].values
        self.g = self.impurity_properties["asm_prm"].values

        assert len(self.mac) == 480 and len(self.ssa) == 480 and len(self.g) == 480


class Ice:
    """Snow or ice column physical properties.

    Instances of Ice contain all the physical properties of each vertical layer of the
    snow or ice column and the underlying surface.

    Attributes:
        dz: array containing thickness of each layer in m
        layer_type: array containing type (0 = grains, 1 = solid ice) in each layer
        cdom: array containing Boolean (1/0) toggling presence of cdom in each layer
        rho: array containing density of each layer in kg/m3
        sfc: array with reflectance of underlying surface per wavelength
        rf: refractive index to use, 0,1,2 or 3 (see docs for definition)
        shp: grain shape per layer where layer_type==1
        rds: grain radius (layer_type==0) or bubble radius (layer_type==0) in each layer
        water: radius of grain+water coating in each layer where layer_type==0
        hex_side: length of each side of hexagonal face for grain_shp==4
        hex_length: column length for hexagonal face for grain_shp==4
        shp_fctr: ratio of nonspherical eff radii to equal vol sphere, in each layer
        ar: aspect ratio of grains in each layer where layer_type==0
        nbr_lyr: number of vertical layers
    """

    def __init__(self, input_file):

        # use config to calculate refractive indices
        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.dz = inputs["ICE"]["DZ"]
        self.layer_type = inputs["ICE"]["LAYER_TYPE"]
        self.cdom = inputs["ICE"]["CDOM"]
        self.rho = inputs["ICE"]["RHO"]
        self.sfc = np.genfromtxt(
            os.getcwd() + "/" + inputs["PATHS"]["SFC"], delimiter="csv"
        )
        self.rf = inputs["ICE"]["RF"]
        self.shp = inputs["ICE"]["SHP"]
        self.rds = inputs["ICE"]["RDS"]
        self.water = inputs["ICE"]["WATER"]
        self.hex_side = inputs["ICE"]["HEX_SIDE"]
        self.hex_length = inputs["ICE"]["HEX_LENGTH"]
        self.shp_fctr = inputs["ICE"]["SHP_FCTR"]
        self.ar = inputs["ICE"]["AR"]
        self.nbr_lyr = len(self.dz)

        self.calculate_refractive_index(input_file)

    def calculate_refractive_index(self, input_file):
        """Calculates ice refractive index from initialized class attributes.

        Takes self.rf and config from inpouts.yaml and uses them to calculate
        new attributes related to the ice refractive index.

        Args:
            self

        Returns:
            ref_idx_im: imaginary part of refractive index
            ref_idx_re: real part of refractive index
            fl_r_dif_a: precomputed diffuse reflectance "perpendicular polarized)
            fl_r_dif_b: precomputed diffuse reflectance "parallel polarized)
            op_dir: directory containing optical properties

        Raises:
            ValueError if rf out of range
        """
        if self.rf < 0 or self.rf > 2:
            raise ValueError("Ice ref index type out of range - between 0 and 2 only")

        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        refidx_file = xr.open_dataset(inputs["PATHS"]["RI_ICE"] + "rfidx_ice.nc")
        fresnel_diffuse_file = xr.open_dataset(
            inputs["PATHS"]["RI_ICE"] + "fl_reflection_diffuse.nc"
        )

        rf = inputs["ICE"]["RF"]
        op_dir_stub = inputs["PATHS"]["OP_DIR_STUBS"][rf]
        ref_idx_name = op_dir_stub[4:9]

        self.ref_idx_re = refidx_file[str("re_" + ref_idx_name)].values
        self.ref_idx_im = refidx_file[str("im_" + ref_idx_name)].values
        self.fl_r_dif_a = fresnel_diffuse_file[
            str("R_dif_fa_ice_" + ref_idx_name)
        ].values
        self.fl_r_dif_b = fresnel_diffuse_file[
            str("R_dif_fb_ice_" + ref_idx_name)
        ].values
        self.op_dir = op_dir_stub


class Illumination:
    """Properties of incoming irradiance.

    Instances of Illumination contain all data relating to the incoming irradiance.

    Attributes:
        direct: Boolean toggling between direct and diffuse irradiance
        solzen: solar zenith angle in degrees from the vertical
        incoming: choice of spectral distribution from file 0-6
        flx_dir: directory containing irradiance files
        stubs: array of stub strings for selecting irradiance files
        nbr_wvl: number fo wavelengths (default 480)
    """

    def __init__(self, input_file):

        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.direct = inputs["RTM"]["DIRECT"]
        self.solzen = inputs["RTM"]["SOLZEN"]
        self.incoming = inputs["RTM"]["INCOMING"]
        self.flx_dir = os.getcwd() + "/" + inputs["PATHS"]["FLX_DIR"]
        self.stubs = inputs["RTM"]["ILLUMINATION_FILE_STUBS"]
        self.nbr_wvl = inputs["RTM"]["NBR_WVL"]

        self.calculate_irradiance()

    def calculate_irradiance(self):
        """Calculates irradiance from initialized attributes.

        Takes mu_not, incoming and file stubs from self and calculates irradiance.

        Args:
            self

        Returns:
            flx_slr: incoming flux from file
            Fd: diffuse irradiance
            Fs: direct irradiance


        Raises:
            ValueError is incoming is out of range
        """

        if self.incoming < 0 or self.incoming > 6:
            raise ValueError("Irradiance type out of range - between 0 and 6 only")

        # update mu_not from solzen
        self.mu_not = np.cos(math.radians(np.rint(self.solzen)))

        # calculate irradiance from file
        cloud_stub = "_cld"
        coszen_stub = ""

        if self.direct:
            cloud_stub = "_clr_"
            coszen_stub = str("SZA" + str(self.solzen).rjust(2, "0"))

        incoming_file = xr.open_dataset(
            str(
                self.flx_dir
                + self.stubs[self.incoming]
                + cloud_stub
                + coszen_stub
                + ".nc"
            )
        )

        flx_slr = incoming_file["flx_frc_sfc"].values
        flx_slr[flx_slr <= 0] = 1e-30
        self.flx_slr = flx_slr
        out = flx_slr / (self.mu_not * np.pi)

        if self.direct:
            self.Fs = out
            self.Fd = np.zeros(self.nbr_wvl)
        else:
            self.Fd = out
            self.Fs = np.zeros(self.nbr_wvl)
        return

        return


class RTConfig:
    """Radiative transfer configuration.

    Attributes:
        aprx_type: choice of two-stream approximation (0-2)
        delta: Boolean to toggle delta transformation (0/1)

    """

    def __init__(self, input_file):

        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.aprx_typ = inputs["RTM"]["APRX_TYP"]
        self.delta = inputs["RTM"]["DELTA"]


class ModelConfig:
    """Model configuration.

    Attributes:
        smooth: Boolean to toggle savitsky-golay filter to smooth albedo
        window_size: window size to use for smoothing func
        poly_order: order of polynomial used to smooth albedo
        dir_base: base directory
        dir_wvl: path to wavelengths in csv file
        sphere_ice_path: directory containing OPs for spherical ice grains
        hex_ice_path: directory containing OPs for hexagonal ice grains
        bubbly_ice_path: directory containing OPs for bubbly ice
        ri_ice_path: path to file containing pure ice refractive index
        op_dir_stubs: sstring stubs for ice optical property files
        wavelengths: array of wavelengths in nm (default 0.205 - 4.995 um)
        nbr_wvl: number of wavelengths (default 480)
        vis_max_idx: index for upper visible wavelength (default 0.75 um)
        nir_max_idx: index for upper NIR wavelength (default 4.995 um)

    """

    def __init__(self, input_file):

        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.smooth = inputs["CTRL"]["SMOOTH"]
        self.window_size = inputs["CTRL"]["WINDOW_SIZE"]
        self.poly_order = inputs["CTRL"]["POLY_ORDER"]
        self.dir_base = str(os.getcwd() + "/")
        self.dir_wvl = inputs["PATHS"]["WVL"]
        self.sphere_ice_path = inputs["PATHS"]["SPHERE_ICE"]
        self.hex_ice_path = inputs["PATHS"]["HEX_ICE"]
        self.bubbly_ice_path = inputs["PATHS"]["BUBBLY_ICE"]
        self.ri_ice_path = inputs["PATHS"]["RI_ICE"]
        self.op_dir_stubs = inputs["PATHS"]["OP_DIR_STUBS"]
        self.wavelengths = np.arange(0.205, 4.999, 0.01)
        self.nbr_wvl = len(self.wavelengths)
        self.vis_max_idx = inputs["RTM"]["VIS_MAX_IDX"]
        self.nir_max_idx = inputs["RTM"]["NIR_MAX_IDX"]


class Outputs:
    """output data from radiative transfer calculations.

    Attributes:
        heat_rt: heating rate in each layer
        BBAVIS: broadband albedo in visible range
        BBANIR: broadband albedo in NIR range
        BBA: broadband albedo across solar spectrum
        abs_slr_btm: absorbed solar energy at bottom surface
        abs_vis_btm: absorbed visible energy at bottom surface
        abs_nir_btm: absorbed NIR energy at bottom surface
        albedo: albedo of ice column
        total_insolation: energy arriving from atmosphere
        abs_slr_tot: total absorbed energy across solar spectrum
        abs_vis_tot: total absorbed energy across visible spectrum
        abs_nir_tot: total absorbed energy across NIR spectrum
        absorbed_flux_per_layer: total absorbed flux per layer
    """

    def __init__(self):
        self.heat_rt = None
        self.BBAVIS = None
        self.BBANIR = None
        self.BBA = None
        self.abs_slr_btm = None
        self.abs_vis_btm = None
        self.abs_nir_btm = None
        self.albedo = None
        self.total_insolation = None
        self.abs_slr_tot = None
        self.abs_vis_tot = None
        self.abs_nir_tot = None
        self.absorbed_flux_per_layer = None


class PlotConfig:
    """Configuration for plotting figures.

    Attributes:
        figsize: size of figure
        facecolor: colour of background
        grid: toggles grid visibility
        grid_color: color of grid lines
        xtick_width: frequency of xticks
        xtick_size: size of ticks on x axis
        ytick_width: frequency of yticks
        ytick_size: size of ticks on y axis
        linewidth: pixel width of line on plot
        fontsize: size of text labels
        xtick_btm: toggles tick positions
        ytick_left: toggle ytick position
        show: toggles showing plot on screen
        save: toggles saving figure to file

    """

    def __init__(self, input_file):

        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.figsize = inputs["PLOT"]["FIG_SIZE"]
        self.facecolor = inputs["PLOT"]["FACECOLOR"]
        self.grid = inputs["PLOT"]["GRID"]
        self.grid_color = inputs["PLOT"]["GRIDCOLOR"]
        self.xtick_width = inputs["PLOT"]["XTICK_WIDTH"]
        self.xtick_size = inputs["PLOT"]["XTICK_SIZE"]
        self.ytick_width = inputs["PLOT"]["YTICK_WIDTH"]
        self.ytick_size = inputs["PLOT"]["YTICK_SIZE"]
        self.linewidth = inputs["PLOT"]["LINEWIDTH"]
        self.fontsize = inputs["PLOT"]["FONTSIZE"]
        self.xtick_btm = inputs["PLOT"]["XTICK_BTM"]
        self.ytick_left = inputs["PLOT"]["YTICK_LEFT"]
        self.save = inputs["PLOT"]["SAVE"]




class BioOpticalConfig:
    
    """Configuration class for bio-optical model.

    Attributes:
        wvl: (numpy array, default: np.arange(0.200, 4.999, 0.001))
                wavelengths in spectral range of interest (in µm, 1nm step)
        wet_density:  (int - used if biomass: True) density of wet biomass
                        (kg/m3 - 1060 and 1160 for snow and glacier algae,
                        Chevrollier et al. 2022)
        dry_density:  (int - used if biomass: True) density of dry biomass
                        (kg/m3 - 625 and 684 for snow and glacier algae,
                        Chevrollier et al. 2022)
        ABS_CFF_CALC: toggles calculating abs_cff from pigments or loading from file.
        abs_cff_loaded_reconstructed: (boolean) True if the
                                abs_cff is loaded as a reconstructed spectrum
                                from pigment absorbance (see methods in
                                Chevrollier et al. 2022)
        abs_cff_loaded_invivo: (boolean) True if the abs_cff is loaded as in vivo
                                spectra of whole cells
        abs_cff_file: (string) directory to the abs_cff file if loaded
        pigment_data: dictionary with pigment file names and associated 
                      intracellular concentrations (ng/cell, ng/µm3 or ng/ng)
        pigment_dir: (string) used if abs_cff_calculated is True, directory to
                         folder containing pigment mass absorption coefficients
                         that must be csv file with size and resolution of wvl,
                         and units in m2/mg
        packaging_correction_SA: (boolean - applied ONLY if
                                  abs_cff_loaded_reconstructed is True) if True,
                                   reconstructed SA abs_cff is corrected for pigment
                                 packaging following Chevrollier et al. 2022
        packaging_correction_GA: (boolean - applied ONLY if
                                abs_cff_loaded_reconstructed is True) if True,
                                reconstructed GA abs_cff is corrected for pigment
                                packaging following Chevrollier et al. 2022
        dir_pckg: (string) directory to pigment packaging correction files
        k_water_dir: (string) path to file with imaginary part of the refractive
                        index of water        
        unit: unit for absorption cross section: 0 = m2/cell, 1 = m2/um3, 3 = m2/mg
              and/or pigment data: 0 = ng/cell, 1 = ng/um3, 3 = ng/mg
        cell_vol: (int - used if cellular: True) volume of the algae
                        cell (um3)
        n_algae: (int) real part of cellular refractive index
                        in the spectral range of wvl (constant 1.38 by default,
                        Chevrollier et al. 2022)
        GO: (boolean) if True, uses geometric optics equations (Cook et
            al. 2020 adapted from Diedenhoven et al (2014)) to calculate single
                scattering OPs assuming cell shape: cylinder
        Mie: (boolean) if True, uses Mie theory to calculate single
                scattering OPs assuming cell shape: sphere
        radius: (int) radius of sphere (Mie)/cynlinder (GO) representing cell (µm)
        length: (int) depth of the cylinder representing the cell (GO option, µm)
        report_dims: (boolean) if True, cell dimensions printed to console
        plot_ssps: (boolean) if True, print plots with ssps
        savefig_ssps: if True, ssps plots saved in the directory
                                savepath
        plot_n_k_abs_cff: (boolean) if True, plot with n,k and abs_cff printed
        saveplots_n_k_abs_cff: (boolean) if True, plots saved in the directory savepath
        savefiles_n_k_abs_cff: (boolean) if True, files with k,n and abs_cff
                        saved in the directory savepath
        savepath: (boolean) directory for saving data if
            savefiles or saveplots toggled on
        smooth: (boolean) if True,  apply optional smoothing filter
        window_size: (int) size of window of smoothing filter (dflt value: 25)
        poly_order: (int) polynomial order of smoothing filter (dflt value: 3)
        save_netcdf = (boolean) if True, saves data in a netcdf file
        savepath_netcdf = (string) save path directory
        filename_netcdf = (string) name of the file containing the
                        optical properties
        information = (string) containing any additional info for
                    metadata in netcdf (e.g. 'Glacier algae OPs derived
                    from GO calculations with empirical abs_cff')

    """

    def __init__(self, input_file):


        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.wvl=np.arange(0.200, 5, 0.001)
        self.wet_density= inputs["BIOOPTICAL"]["WET_DENSITY"]
        self.dry_density= inputs["BIOOPTICAL"]["DRY_DENSITY"]
        self.abs_cff_calculated=inputs["BIOOPTICAL"]["ABS_CFF_CALC"]
        self.abs_cff_loaded_reconstructed=inputs["BIOOPTICAL"]["ABS_CFF_LOAD_RECONSTRUCTED"]
        self.abs_cff_loaded_invivo=inputs["BIOOPTICAL"]["ABS_CFF_LOAD_INVIVO"]
        self.abs_cff_file=inputs["BIOOPTICAL"]["ABS_CFF_FILE"]
        self.pigment_data=inputs["BIOOPTICAL"]["PIGMENT_CONC"]
        self.pigment_dir=inputs["BIOOPTICAL"]["PIGMENT_DIR"]
        self.packaging_correction_SA=inputs["BIOOPTICAL"]["PCKG_SA"]
        self.packaging_correction_GA=inputs["BIOOPTICAL"]["PCKG_GA"]
        self.dir_pckg=inputs["BIOOPTICAL"]["DIR_PCKG"]
        self.k_water_dir=inputs["BIOOPTICAL"]["K_WATER_DIR"]
        self.unit=inputs["BIOOPTICAL"]["UNIT"]
        self.cell_vol=inputs["BIOOPTICAL"]["CELL_VOLUME"]
        self.n_algae=inputs["BIOOPTICAL"]["N_ALGAE"]
        self.GO = inputs["BIOOPTICAL"]["GO"]
        self.Mie = inputs["BIOOPTICAL"]["MIE"]
        self.radius = inputs["BIOOPTICAL"]["CELL_R"]
        self.length = inputs["BIOOPTICAL"]["CELL_L"]
        self.report_dims = inputs["BIOOPTICAL"]["REPORT_DIMS"] 
        self.plot_ssps = inputs["BIOOPTICAL"]["PLOT_SSPS"]
        self.savefig_ssps = inputs["BIOOPTICAL"]["SAVEFIG_SSPS"]
        self.plot_k_abs_cff = inputs["BIOOPTICAL"]["PLOT_K_ABS_CFF"]
        self.saveplots_k_abs_cff = inputs["BIOOPTICAL"]["SAVE_PLOT_K_ABS_CFF"]
        self.savefiles_n_k_abs_cff = inputs["BIOOPTICAL"]["SAVE_N_K_ABS_CFF"]
        self.savepath = inputs["BIOOPTICAL"]["SAVE_PATH_FIG"]
        self.smooth = inputs["CTRL"]["SMOOTH"]
        self.window_size = inputs["CTRL"]["WINDOW_SIZE"]
        self.poly_order = inputs["CTRL"]["POLY_ORDER"]
        self.save_netcdf = inputs["BIOOPTICAL"]["SAVE_NETCDF"]
        self.savepath_netcdf = inputs["BIOOPTICAL"]["SAVE_PATH_NETCDF"]
        self.filename_netcdf = inputs["BIOOPTICAL"]["FILENAME_NETCDF"]
        self.information = inputs["BIOOPTICAL"]["INFO_NETCDF"]

        self.validate_biooptical_inputs()

    def validate_biooptical_inputs(self):
        
        if self.Mie:
            assert(self.GO==False) 
        if self.GO:
            assert(self.Mie==False)  
          
        return


if __name__ == '__main__':
    pass