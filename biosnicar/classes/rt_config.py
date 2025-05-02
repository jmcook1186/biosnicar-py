import yaml

class RTConfig:
    """Radiative transfer configuration.

    Attributes:
        aprx_type: choice of two-stream approximation (0-2)
        delta: Boolean to toggle delta transformation (0/1)

    """

    def __init__(self, input_file):
        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.aprx_typ = inputs["RTM"]["APRX_TYP"]
        self.delta = inputs["RTM"]["DELTA"] 