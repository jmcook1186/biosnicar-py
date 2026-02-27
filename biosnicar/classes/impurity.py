import os
import numpy as np
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
        file: name of file containing optical properties and size distribution
        mac: mass absorption coefficient (m2/kg or m2/cell)
        ssa: single scattering albedo
        g: asymmetry parameter

    """

    def __init__(self, file, coated, unit, name, conc):
        self.name = name
        self.unit = unit
        self.conc = conc
        self.file = file

        npz_file = os.path.splitext(file)[0] + ".npz"
        npz_path = os.path.join(
            os.path.dirname(os.path.dirname(biosnicar.__file__)),
            "Data", "OP_data", "480band", "lap", npz_file,
        )
        impurity_properties = np.load(npz_path)

        if coated:
            mac_stub = "ext_cff_mss_ncl"
        elif (name == "ga") or (name == "sa"):
            mac_stub = "ext_xsc"
        else:
            mac_stub = "ext_cff_mss"

        self.mac = impurity_properties[mac_stub]
        self.ssa = impurity_properties["ss_alb"]
        self.g = impurity_properties["asm_prm"]

        assert len(self.mac) == 480 and len(self.ssa) == 480 and len(self.g) == 480