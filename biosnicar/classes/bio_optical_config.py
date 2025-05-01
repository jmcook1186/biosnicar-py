import numpy as np
import yaml

class BioOpticalConfig:
    """Configuration class for bio-optical model.

    Attributes:
        wvl: (numpy array, default: np.arange(0.200, 4.999, 0.001))
            wavelengths in spectral range of interest (in µm, 1nm step)
        wet_density:  (int - used if biomass: True) density of wet biomass
        (kg/m3 - 1060 and 1160 for snow and glacier algae,Chevrollier et al. 2022)
        dry_density:  (int - used if biomass: True) density of dry biomass
            (kg/m3 - 625 and 684 for snow and glacier algae,
            Chevrollier et al. 2022)
        ABS_CFF_CALC: toggles calculating abs_cff from pigments or loading from file.
        abs_cff_loaded_reconstructed: (boolean) True if the
            abs_cff is loaded as a reconstructed spectrum
            from pigment absorbance (see methods in Chevrollier et al. 2022)
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
        cell_vol: (int - used if cellular: True) volume of the algae cell (um3)
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
        savefig_ssps: if True, ssps plots saved in the directory savepath
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

        self.wvl = np.arange(0.200, 5, 0.001)
        self.wet_density = inputs["BIOOPTICAL"]["WET_DENSITY"]
        self.dry_density = inputs["BIOOPTICAL"]["DRY_DENSITY"]
        self.abs_cff_calculated = inputs["BIOOPTICAL"]["ABS_CFF_CALC"]
        self.abs_cff_loaded_reconstructed = inputs["BIOOPTICAL"][
            "ABS_CFF_LOAD_RECONSTRUCTED"
        ]
        self.abs_cff_loaded_invivo = inputs["BIOOPTICAL"]["ABS_CFF_LOAD_INVIVO"]
        self.abs_cff_file = inputs["BIOOPTICAL"]["ABS_CFF_FILE"]
        self.pigment_data = inputs["BIOOPTICAL"]["PIGMENT_CONC"]
        self.pigment_dir = inputs["BIOOPTICAL"]["PIGMENT_DIR"]
        self.packaging_correction_SA = inputs["BIOOPTICAL"]["PCKG_SA"]
        self.packaging_correction_GA = inputs["BIOOPTICAL"]["PCKG_GA"]
        self.dir_pckg = inputs["BIOOPTICAL"]["DIR_PCKG"]
        self.k_water_dir = inputs["BIOOPTICAL"]["K_WATER_DIR"]
        self.unit = inputs["BIOOPTICAL"]["UNIT"]
        self.cell_vol = inputs["BIOOPTICAL"]["CELL_VOLUME"]
        self.n_algae = inputs["BIOOPTICAL"]["N_ALGAE"]
        self.GO = inputs["BIOOPTICAL"]["GO"]
        self.Mie = inputs["BIOOPTICAL"]["MIE"]
        self.radius = inputs["BIOOPTICAL"]["CELL_R"]
        self.length = inputs["BIOOPTICAL"]["CELL_L"]
        self.report_dims = inputs["BIOOPTICAL"]["REPORT_DIMS"]
        self.plot_ssps = inputs["BIOOPTICAL"]["PLOT_SSPS"]
        self.savefig_ssps = inputs["BIOOPTICAL"]["SAVEFIG_SSPS"]
        self.plot_n_k_abs_cff = inputs["BIOOPTICAL"]["PLOT_N_K_ABS_CFF"]
        self.saveplots_n_k_abs_cff = inputs["BIOOPTICAL"]["SAVEPLOTS_N_K_ABS_CFF"]
        self.savefiles_n_k_abs_cff = inputs["BIOOPTICAL"]["SAVEFILES_N_K_ABS_CFF"]
        self.savepath = inputs["BIOOPTICAL"]["SAVEPATH"]
        self.smooth = inputs["BIOOPTICAL"]["SMOOTH"]
        self.window_size = inputs["BIOOPTICAL"]["WINDOW_SIZE"]
        self.poly_order = inputs["BIOOPTICAL"]["POLY_ORDER"]
        self.save_netcdf = inputs["BIOOPTICAL"]["SAVE_NETCDF"]
        self.savepath_netcdf = inputs["BIOOPTICAL"]["SAVEPATH_NETCDF"]
        self.filename_netcdf = inputs["BIOOPTICAL"]["FILENAME_NETCDF"]
        self.information = inputs["BIOOPTICAL"]["INFO_NETCDF"] 