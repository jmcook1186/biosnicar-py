"""Optical properties calculation modules for BioSNICAR.

This package contains modules for calculating various optical properties:
- Bubble effective radius calculation
- Column optical properties
- Geometric optics for ice
- Mie theory for coated water spheres
- Single scattering properties for spheres
"""

from biosnicar.optical_properties.bubble_reff_calculator import *
from biosnicar.optical_properties.column_OPs import *
from biosnicar.optical_properties.geometric_optics_ice import *
from biosnicar.optical_properties.mie_coated_water_spheres import *
from biosnicar.optical_properties.ssps_spheres_generator import *

__all__ = [
    'bubble_reff_calculator',
    'column_OPs',
    'geometric_optics_ice',
    'mie_coated_water_spheres',
    'ssps_spheres_generator'
] 