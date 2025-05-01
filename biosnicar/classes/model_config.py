import os
import numpy as np
import yaml
import biosnicar

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
        self.dir_base = str(os.path.dirname(os.path.dirname(biosnicar.__file__)))+ "/"
        self.dir_wvl = inputs["PATHS"]["WVL"]
        self.sphere_ice_path = inputs["PATHS"]["SPHERE_ICE"]
        self.sphere_water_path = inputs["PATHS"]["SPHERE_WATER"]
        self.fn_ice = inputs["PATHS"]["FN_ICE"]
        self.fn_water = inputs["PATHS"]["FN_WATER"]
        self.hex_ice_path = inputs["PATHS"]["HEX_ICE"]
        self.bubbly_ice_path = inputs["PATHS"]["BUBBLY_ICE"]
        self.ri_ice_path = inputs["PATHS"]["RI_ICE"]
        self.op_dir_stubs = inputs["PATHS"]["OP_DIR_STUBS"]
        self.savefigpath = inputs["PLOT"]["SAVEPATH"]
        self.wavelengths = np.arange(0.205, 4.999, 0.01)
        self.nbr_wvl = len(self.wavelengths)
        self.vis_max_idx = inputs["RTM"]["VIS_MAX_IDX"]
        self.nir_max_idx = inputs["RTM"]["NIR_MAX_IDX"] 