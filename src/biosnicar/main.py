#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from pathlib import Path
from get_albedo import get_albedo

# call easy albedo func
albedo = get_albedo("toon", plot=True, validate=True)
