#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pathlib import Path
from biosnicar import get_albedo

# call easy albedo func
albedo = get_albedo.get("adding-doubling", plot=True, validate=True)
