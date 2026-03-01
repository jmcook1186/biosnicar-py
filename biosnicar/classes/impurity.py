import os
import numpy as np
import biosnicar

# Module-level cache: loaded once on first use
_lap_data = None


def _get_lap_data():
    global _lap_data
    if _lap_data is None:
        npz_path = str(biosnicar.DATA_DIR / "OP_data" / "480band" / "lap.npz")
        _lap_data = np.load(npz_path)
    return _lap_data


def invalidate_lap_cache():
    """Reset the module-level LAP cache so new data is picked up."""
    global _lap_data
    _lap_data = None


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

        stem = os.path.splitext(file)[0]
        lap = _get_lap_data()

        if coated:
            mac_stub = "ext_cff_mss_ncl"
        elif (name == "ga") or (name == "sa"):
            mac_stub = "ext_xsc"
        else:
            mac_stub = "ext_cff_mss"

        self.mac = lap[f"{stem}__{mac_stub}"]
        self.ssa = lap[f"{stem}__ss_alb"]
        self.g = lap[f"{stem}__asm_prm"]

        assert len(self.mac) == 480 and len(self.ssa) == 480 and len(self.g) == 480
