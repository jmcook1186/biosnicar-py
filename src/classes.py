from pathlib import Path
import xarray as xr
import yaml
import numpy as np
import math


class Impurity:
    def __init__(self, dir_base, file, coated, cfactor, unit, name, conc):

        self.name = name
        self.cfactor = cfactor
        self.unit = unit
        self.conc = conc
        self.file = file

        self.impurity_properties = xr.open_dataset(
            str(dir_base + "/Data/OP_data/480band/lap/" + file)
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
    def __init__(self):

        # use config to calculate refractive indices
        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.dz = inputs["VARIABLES"]["DZ"]
        self.layer_type = inputs["VARIABLES"]["LAYER_TYPE"]
        self.cdom = inputs["VARIABLES"]["CDOM"]
        self.rho = inputs["VARIABLES"]["RHO"]
        self.sfc = np.genfromtxt(
            inputs["PATHS"]["DIR_BASE"] + inputs["PATHS"]["SFC"], delimiter="csv"
        )
        self.rf = inputs["VARIABLES"]["RF"]
        self.shp = inputs["VARIABLES"]["SHP"]
        self.rds = inputs["VARIABLES"]["RDS"]
        self.water = inputs["VARIABLES"]["WATER"]
        self.hex_side = inputs["VARIABLES"]["HEX_SIDE"]
        self.hex_length = inputs["VARIABLES"]["HEX_LENGTH"]
        self.shp_fctr = inputs["VARIABLES"]["SHP_FCTR"]
        self.ar = inputs["VARIABLES"]["AR"]
        self.nbr_lyr = len(self.dz)

        self.calculate_refractive_index()

    def calculate_refractive_index(self):

        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        refidx_file = xr.open_dataset(inputs["PATHS"]["RI_ICE"] + "rfidx_ice.nc")
        fresnel_diffuse_file = xr.open_dataset(
            inputs["PATHS"]["RI_ICE"] + "fl_reflection_diffuse.nc"
        )

        rf = inputs["VARIABLES"]["RF"]
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
    def __init__(self):

        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.direct = inputs["RTM"]["DIRECT"]
        self.solzen = inputs["RTM"]["SOLZEN"]
        self.incoming = inputs["RTM"]["INCOMING"]
        self.flx_dir = inputs["PATHS"]["DIR_BASE"] + inputs["PATHS"]["FLX_DIR"]
        self.stubs = inputs["RTM"]["ILLUMINATION_FILE_STUBS"]
        self.nbr_wvl = inputs["RTM"]["NBR_WVL"]

        self.calculate_irradiance()
        
    def calculate_irradiance(self):

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
    def __init__(self):

        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.aprx_typ = inputs["RTM"]["APRX_TYP"]
        self.delta = inputs["RTM"]["DELTA"]



class ModelConfig:
    def __init__(self):

        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.smooth = inputs["CTRL"]["SMOOTH"]
        self.window_size = inputs["CTRL"]["WINDOW_SIZE"]
        self.poly_order = inputs["CTRL"]["POLY_ORDER"]
        self.smooth = inputs["CTRL"]["SMOOTH"]
        self.dir_base = inputs["PATHS"]["DIR_BASE"]
        self.dir_wvl = inputs["PATHS"]["WVL"]
        self.sphere_ice_path = inputs["PATHS"]["SPHERE_ICE"]
        self.hex_ice_path = inputs["PATHS"]["HEX_ICE"]
        self.bubbly_ice_path = inputs["PATHS"]["BUBBLY_ICE"]
        self.ri_ice_path = inputs["PATHS"]["RI_ICE"]
        self.sphere_ice_path = inputs["PATHS"]["SPHERE_ICE"]
        self.sphere_ice_path = inputs["PATHS"]["SPHERE_ICE"]
        self.op_dir_stubs = inputs["PATHS"]["OP_DIR_STUBS"]
        self.wavelengths = np.arange(0.205, 4.999, 0.01)
        self.nbr_wvl = len(self.wavelengths)
        self.vis_max_idx = inputs["RTM"]["VIS_MAX_IDX"]
        self.nir_max_idx = inputs["RTM"]["NIR_MAX_IDX"]


class Outputs:
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
    def __init__(self):

        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # self.fig_width = inputs["PLOT"]["FIG_WIDTH"]
        # self.fig_height = inputs["PLOT"]["FIG_HEIGHT"]
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
        self.show = inputs["PLOT"]["SHOW"]
        self.save = inputs["PLOT"]["SAVE"]


class DisplayConfig:
    def __init__(self):
        with open("./src/inputs.yaml", "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

            self.print_bba = inputs["CTRL"]["PRINT_BBA"]
            self.print_band_ratios = inputs["CTRL"]["PRINT_BAND_RATIOS"]


