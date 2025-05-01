import os
import numpy as np
import pandas as pd
import xarray as xr
import yaml
import biosnicar

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
            str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/" + inputs["PATHS"]["SFC"]), delimiter="csv"
        )
        self.rf = inputs["ICE"]["RF"]
        self.shp = inputs["ICE"]["SHP"]
        self.rds = inputs["ICE"]["RDS"]
        self.water = inputs["ICE"]["WATER_COATING"]
        self.hex_side = inputs["ICE"]["HEX_SIDE"]
        self.hex_length = inputs["ICE"]["HEX_LENGTH"]
        self.shp_fctr = inputs["ICE"]["SHP_FCTR"]
        self.ar = inputs["ICE"]["AR"]
        self.nbr_lyr = len(self.dz)
        self.lwc = inputs["ICE"]["LWC"]
        self.lwc_pct_bbl = inputs["ICE"]["LWC_PCT_BBL"]
        self.ref_idx_im_water = pd.read_csv(
            str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/" + inputs["PATHS"]["RI_ICE"] + 
            'refractive_index_water_273K_Rowe2020.csv'
        )).k.values
    
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

        refidx_file = xr.open_dataset(
            str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/" + inputs["PATHS"]["RI_ICE"] + "rfidx_ice.nc"))
        fresnel_diffuse_file = xr.open_dataset(
            str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/" +
            inputs["PATHS"]["RI_ICE"] + "fl_reflection_diffuse.nc"))

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