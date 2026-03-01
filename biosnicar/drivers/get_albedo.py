#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from biosnicar.drivers.run_model import run_model


def get(solver, plot, validate):
    """Run BioSNICAR with default inputs and return spectral albedo.

    Thin wrapper around :func:`~biosnicar.drivers.run_model.run_model`
    that returns only the albedo array for backwards compatibility.

    Args:
        solver: ``"adding-doubling"`` or ``"toon"``.
        plot: If True, display a spectral albedo plot.
        validate: If True, validate inputs before running.

    Returns:
        numpy array of spectral albedo (480 wavelengths).
    """
    outputs = run_model(solver=solver, validate=validate, plot=plot)
    return outputs.albedo
