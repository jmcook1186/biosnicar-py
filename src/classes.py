from pathlib import Path
import xarray as xr
import yaml
import numpy as np
import math

dir_base = "/home/joe/Code/BioSNICAR_GO_PY/Data/OP_data/480band/"

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
    def __init__(self, dir_base):
        with open("/home/joe/Code/BioSNICAR_GO_PY/src/ice_physical_config.yaml", "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
        self.dz = cfg["dz"]
        self.layer_type = cfg["layer_type"]
        self.cdom = cfg["cdom"]
        self.rho = cfg["rho"]        
        self.sfc = np.genfromtxt(dir_base+cfg["sfc_file"],delimiter="csv")
        self.rf = cfg["rf"]
        self.shp = cfg["shp"]
        self.rds = cfg["rds"]
        self.water = cfg["water"]
        self.hex_side = cfg["hex_side"]
        self.hex_length = cfg["hex_length"]
        self.shp_fctr = cfg["shp_fctr"]
        self.ar = cfg["ar"]



class illumination:
    def __init__(self):
    
        with open("/home/joe/Code/BioSNICAR_GO_PY/src/rtm_config.yaml", "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
        self.direct = cfg["direct"]
        self.solzen = cfg["solzen"]
        self.incoming = cfg["incoming"]
        self.mu_not = np.cos(math.radians(np.rint(self.solzen)))


    def calculate_flx_slr(self):

        with open("/home/joe/Code/BioSNICAR_GO_PY/src/rtm_config.yaml", "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        
            stubs = cfg["illumination_file_stubs"]
            flx_dir = cfg["flx_dir"]
            nbr_wvl = cfg["nbr_wvl"]
        
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
    def __init__(self):

        with open("/home/joe/Code/BioSNICAR_GO_PY/src/rtm_config.yaml", "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.aprx_typ = cfg["aprx_typ"]
        self.delta = cfg["delta"]
        self.toon = cfg["toon"]
        self.add_double = cfg["add_double"]


