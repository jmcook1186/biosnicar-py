from biosnicar.utils.load_inputs import load_inputs

class RTConfig:
    """Radiative transfer configuration.

    Attributes:
        aprx_type: choice of two-stream approximation (0-2)
        delta: Boolean to toggle delta transformation (0/1)

    """

    def __init__(self, input_file):
        inputs = load_inputs(input_file)

        self.aprx_typ = inputs["RTM"]["APRX_TYP"]
        self.delta = inputs["RTM"]["DELTA"] 