#!/usr/bin/python
"""Calculates effective radius of air bubbles.

This script is for taking specific surface area for bubbly ice and calculating the
effective radius of the bubbles assuming a lognormal bubble size distribution.

The reason this is required is that the well-known conversion between SSA and r_eff:

r_eff = 3/(P_i * SSA) where P_i is density of ice (917 kg m-3)

gives the effective radius of a discrete grain of given SSA, for a collection of bubbles
in a bulk medium of ice a more nuanced calculation is required.

The derivation of the calculations are explained in the pdf ./assets/SSA_derivation.pdf
from Chloe Whicker (UMich).
"""
