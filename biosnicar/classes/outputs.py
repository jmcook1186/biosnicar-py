class Outputs:
    """output data from radiative transfer calculations.

    Attributes:
        heat_rt: heating rate in each layer
        BBAVIS: broadband albedo in visible range
        BBANIR: broadband albedo in NIR range
        BBA: broadband albedo across solar spectrum
        abs_slr_btm: absorbed solar energy at bottom surface
        abs_vis_btm: absorbed visible energy at bottom surface
        abs_nir_btm: absorbed NIR energy at bottom surface
        albedo: albedo of ice column
        total_insolation: energy arriving from atmosphere
        abs_slr_tot: total absorbed energy across solar spectrum
        abs_vis_tot: total absorbed energy across visible spectrum
        abs_nir_tot: total absorbed energy across NIR spectrum
        absorbed_flux_per_layer: total absorbed flux per layer
    """

    def __init__(self):
        self.heat_rt = None
        self.BBAVIS = None
        self.BBANIR = None
        self.BBA = None
        self.abs_slr_btm = None
        self.abs_vis_btm = None
        self.abs_nir_btm = None
        self.albedo = None
        self.total_insolation = None
        self.abs_slr_tot = None
        self.abs_vis_tot = None
        self.abs_nir_tot = None
        self.absorbed_flux_per_layer = None 