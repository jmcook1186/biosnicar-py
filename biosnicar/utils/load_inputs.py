"""Cached YAML loader for input configuration files.

Avoids redundant parsing when multiple classes are constructed from the
same inputs file during a single model setup.
"""

import yaml

_cache = {}


def load_inputs(input_file):
    """Load and cache a YAML input file.

    Returns the parsed dict, reusing a cached copy if the same *input_file*
    path has already been loaded.

    Args:
        input_file: path to the YAML configuration file.

    Returns:
        Parsed dict of the YAML contents.
    """
    if input_file not in _cache:
        with open(input_file, "r") as ymlfile:
            _cache[input_file] = yaml.load(ymlfile, Loader=yaml.FullLoader)
    return _cache[input_file]


def clear_cache():
    """Clear the YAML input cache.

    Call this if you need to reload an input file after modifying it on disk.
    """
    _cache.clear()
