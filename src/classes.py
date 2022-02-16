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

        self.impurity_properties = xr.open_dataset(str(dir_base + file))

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
        
        self.dz = ice_cfg["dz"]
        self.layer_type = ice_cfg["layer_type"]
        self.cdom = ice_cfg["cdom"]
        self.rho = ice_cfg["rho"]        
        self.sfc = np.genfromtxt(model_cfg["DIR_BASE"]+ice_cfg["sfc_file"],delimiter="csv")
        self.rf = ice_cfg["rf"]
        self.shp = ice_cfg["shp"]
        self.rds = ice_cfg["rds"]
        self.water = ice_cfg["water"]
        self.hex_side = ice_cfg["hex_side"]
        self.hex_length = ice_cfg["hex_length"]
        self.shp_fctr = ice_cfg["shp_fctr"]
        self.ar = ice_cfg["ar"]



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


