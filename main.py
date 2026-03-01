#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from biosnicar.drivers import get_albedo

# call easy albedo func
albedo = get_albedo.get("adding-doubling", plot=True, validate=True)
