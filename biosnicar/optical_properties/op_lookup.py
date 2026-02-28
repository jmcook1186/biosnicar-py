"""Lookup table loaders for pre-computed optical properties.

Provides O(1) access to optical property arrays stored in compressed
.npz files, replacing per-file netCDF reads.
"""

import numpy as np


# Module-level cache: {npz_path: OpLookupTable or HexOpLookupTable}
_cache = {}


def get_lut(path):
    """Return a cached OpLookupTable for the given .npz path."""
    if path not in _cache:
        _cache[path] = OpLookupTable(path)
    return _cache[path]


def get_hex_lut(path):
    """Return a cached HexOpLookupTable for the given .npz path."""
    if path not in _cache:
        _cache[path] = HexOpLookupTable(path)
    return _cache[path]


class OpLookupTable:
    """Lookup table for 1D-indexed (radius) optical properties.

    Loads a .npz file containing arrays indexed by radius. Provides
    O(1) lookup by radius for any stored variable.

    Attributes:
        data: loaded .npz archive
        _radius_to_idx: dict mapping integer radius to array index
    """

    def __init__(self, path):
        self.data = np.load(str(path))
        radii = self.data["radii"]
        self._radius_to_idx = {int(r): i for i, r in enumerate(radii)}

    def get(self, radius, var_name):
        """Look up a 480-element array by radius and variable name.

        Args:
            radius: integer radius in micrometers
            var_name: variable name (e.g. 'ss_alb', 'ext_cff_mss')

        Returns:
            numpy array of shape (480,)
        """
        idx = self._radius_to_idx[int(radius)]
        return self.data[var_name][idx]


class HexOpLookupTable:
    """Lookup table for 2D-indexed (side, length) optical properties.

    Loads a .npz file containing arrays indexed by (side, length) pairs
    for hexagonal ice column optical properties.

    Attributes:
        data: loaded .npz archive
        _key_to_idx: dict mapping (side, length) tuple to array index
    """

    def __init__(self, path):
        self.data = np.load(str(path))
        sides = self.data["sides"]
        lengths = self.data["lengths"]
        self._key_to_idx = {
            (int(s), int(l)): i for i, (s, l) in enumerate(zip(sides, lengths))
        }

    def get(self, side, length, var_name):
        """Look up a 480-element array by side, length, and variable name.

        Args:
            side: integer side length in micrometers
            length: integer column length in micrometers
            var_name: variable name (e.g. 'ss_alb', 'ext_cff_mss')

        Returns:
            numpy array of shape (480,)
        """
        idx = self._key_to_idx[(int(side), int(length))]
        return self.data[var_name][idx]
