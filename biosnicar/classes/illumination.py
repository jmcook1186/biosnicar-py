import os
import math
import numpy as np
import xarray as xr
import yaml
import biosnicar

class Illumination:
    """Properties of incoming irradiance.

    Instances of Illumination contain all data relating to the incoming irradiance.

    Attributes:
        direct: Boolean toggling between direct and diffuse irradiance
        solzen: solar zenith angle in degrees from the vertical
        incoming: choice of spectral distribution from file 0-6
        flx_dir: directory containing irradiance files
        stubs: array of stub strings for selecting irradiance files
        nbr_wvl: number fo wavelengths (default 480)
    """

    def __init__(self, input_file):
        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.direct = inputs["RTM"]["DIRECT"]
        self.solzen = inputs["RTM"]["SOLZEN"]
        self.incoming = inputs["RTM"]["INCOMING"]
        self.flx_dir = str(os.path.dirname(os.path.dirname(biosnicar.__file__))
                 + "/" + inputs["PATHS"]["FLX_DIR"])
        self.stubs = inputs["PATHS"]["ILLUMINATION_FILE_STUBS"]
        self.nbr_wvl = inputs["RTM"]["NBR_WVL"]

        self.calculate_irradiance()

    def calculate_irradiance(self):
        """Calculates irradiance from initialized attributes.

        Takes mu_not, incoming and file stubs from self and calculates irradiance.

        Args:
            self

        Returns:
            flx_slr: incoming flux from file
            Fd: diffuse irradiance
            Fs: direct irradiance


        Raises:
            ValueError is incoming is out of range
        """

        if self.incoming < 0 or self.incoming > 6:
            raise ValueError("Irradiance type out of range - between 0 and 6 only")

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