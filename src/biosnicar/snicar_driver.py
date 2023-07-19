#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pathlib import Path
from get_albedo import get_albedo


if __name__ == "__main__":
    BIOSNICAR_SRC_PATH = Path(__file__).resolve().parent

    # define input file
    input_file = BIOSNICAR_SRC_PATH.joinpath("inputs.yaml").as_posix()
    out = get_albedo(input_file, "adding-doubling")
