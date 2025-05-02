import os
import numpy as np
import xarray as xr
import biosnicar

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

    def __init__(self, file, coated, unit, name, conc):
        self.name = name
        self.unit = unit
        self.conc = conc
        self.file = file

        self.impurity_properties = xr.open_dataset(
            str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/Data/OP_data/480band/lap/" + file)
        )

        if coated:
            mac_stub = "ext_cff_mss_ncl"
        elif (name == "ga") or (name == "sa"):
            mac_stub = "ext_xsc"
        else:
            mac_stub = "ext_cff_mss"

        self.mac = self.impurity_properties[mac_stub].values
        self.ssa = self.impurity_properties["ss_alb"].values
        self.g = self.impurity_properties["asm_prm"].values

        assert len(self.mac) == 480 and len(self.ssa) == 480 and len(self.g) == 480 