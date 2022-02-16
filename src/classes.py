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

        self.impurity_properties = xr.open_dataset(str(dir_base + "/Data/OP_data/480band/lap/" + file))

        if coated:
            mac_stub = "ext_cff_mss_ncl"
        else:
            mac_stub = "ext_cff_mss"

        self.mac = self.impurity_properties[mac_stub].values
        self.ssa = self.impurity_properties["ss_alb"].values
        self.g = self.impurity_properties["asm_prm"].values

        assert (len(self.mac) == 480 and len(self.ssa) == 480 and len(self.g) == 480)


class Ice:

    def __init__(self, ice_cfg, model_cfg):
        
        self.dz = ice_cfg["VARIABLES"]["dz"]
        self.layer_type = ice_cfg["VARIABLES"]["layer_type"]
        self.cdom = ice_cfg["VARIABLES"]["cdom"]
        self.rho = ice_cfg["VARIABLES"]["rho"]        
        self.sfc = np.genfromtxt(model_cfg["PATHS"]["DIR_BASE"]+ice_cfg["PATHS"]["SFC"],delimiter="csv")
        self.rf = ice_cfg["VARIABLES"]["rf"]
        self.shp = ice_cfg["VARIABLES"]["shp"]
        self.rds = ice_cfg["VARIABLES"]["rds"]
        self.water = ice_cfg["VARIABLES"]["water"]
        self.hex_side = ice_cfg["VARIABLES"]["hex_side"]
        self.hex_length = ice_cfg["VARIABLES"]["hex_length"]
        self.shp_fctr = ice_cfg["VARIABLES"]["shp_fctr"]
        self.ar = ice_cfg["VARIABLES"]["ar"]
        self.nbr_lyr = len(self.dz)


    def get_ref_indices(self, ice_cfg):

        refidx_file = xr.open_dataset(ice_cfg["PATHS"]["RI_ICE"] + "rfidx_ice.nc")
        fresnel_diffuse_file = xr.open_dataset(ice_cfg["PATHS"]["RI_ICE"] + "fl_reflection_diffuse.nc")

        rf = ice_cfg["VARIABLES"]["rf"]
        op_dir_stub = ice_cfg["PATHS"]["OP_DIR_STUBS"][rf]
        ref_idx_name = op_dir_stub[4:9]
        self. ref_idx_re = refidx_file[str("re_" + ref_idx_name)].values
        self. ref_idx_im = refidx_file[str("im_" + ref_idx_name)].values
        self.fl_r_dif_a = fresnel_diffuse_file[str("R_dif_fa_ice_" + ref_idx_name)].values
        self.fl_r_dif_b = fresnel_diffuse_file[str("R_dif_fb_ice_" + ref_idx_name)].values
        self.op_dir = op_dir_stub




class Illumination:
    def __init__(self, rtm_cfg):
    
        self.direct = rtm_cfg["direct"]
        self.solzen = rtm_cfg["solzen"]
        self.incoming = rtm_cfg["incoming"]
        self.mu_not = np.cos(math.radians(np.rint(self.solzen)))
        self.rtm_cfg = rtm_cfg

    def calculate_flx_slr(self):

        stubs = self.rtm_cfg["illumination_file_stubs"]
        flx_dir = self.rtm_cfg["flx_dir"]
        nbr_wvl = self.rtm_cfg["nbr_wvl"]
        
        cloud_stub = "_cld"
        coszen_stub = ""
        
        if self.direct:
            cloud_stub = "_clr_"
            coszen_stub = str("SZA" + str(self.solzen).rjust(2, "0"))

        incoming_file = xr.open_dataset(
            str(flx_dir + stubs[self.incoming] + cloud_stub + coszen_stub + ".nc"))

        flx_slr = incoming_file["flx_dwn_sfc"].values
        flx_slr[flx_slr <= 0] = 1e-30
        flx_slr = flx_slr
        out = flx_slr / (self.mu_not * np.pi)
        
        if self.direct:
            self.Fs = out
            self.Fd = np.zeros(nbr_wvl)
        else:
            self.Fd = out
            self.Fs = np.zeros(nbr_wvl)

        return




class RTConfig:
    def __init__(self, rtm_cfg):

        self.aprx_typ = rtm_cfg["aprx_typ"]
        self.delta = rtm_cfg["delta"]
        self.toon = rtm_cfg["toon"]
        self.add_double = rtm_cfg["add_double"]


