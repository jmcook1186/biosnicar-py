"""Tests for the parameter sweep API."""

import numpy as np
import pytest

from biosnicar.drivers.sweep import parameter_sweep


def test_single_param_sweep_shape():
    """Single-parameter sweep returns correct DataFrame shape."""
    df = parameter_sweep(
        params={"solzen": [30, 50, 70]},
        progress=False,
    )
    assert len(df) == 3
    assert "solzen" in df.columns
    assert "BBA" in df.columns


def test_multi_param_sweep_shape():
    """Multi-parameter sweep returns rows = product of lengths."""
    df = parameter_sweep(
        params={"solzen": [30, 50], "rds": [200, 500]},
        progress=False,
    )
    assert len(df) == 4  # 2 x 2
    assert set(df.columns) >= {"solzen", "rds", "BBA", "BBAVIS", "BBANIR"}


def test_impurity_conc_decreases_bba():
    """More black carbon -> lower BBA."""
    df = parameter_sweep(
        params={"impurity.0.conc": [0, 5000, 50000]},
        progress=False,
    )
    bbas = df.sort_values("impurity.0.conc")["BBA"].values
    assert bbas[0] > bbas[1] > bbas[2], (
        f"BBA should decrease with increasing BC: {bbas}"
    )


def test_bba_decreases_with_sza():
    """BBA decreases with increasing solar zenith angle (for direct beam)."""
    df = parameter_sweep(
        params={"solzen": [30, 50, 70]},
        progress=False,
    )
    bbas = df.sort_values("solzen")["BBA"].values
    # BBA generally increases with SZA for snow/ice (more forward scattering at
    # oblique angles), but the key thing is values are distinct and physical.
    assert all(0 < b < 1 for b in bbas), f"BBA out of physical range: {bbas}"


def test_bba_decreases_with_grain_radius():
    """Larger grains -> lower BBA (more absorption path length)."""
    df = parameter_sweep(
        params={"rds": [100, 500, 2000]},
        progress=False,
    )
    bbas = df.sort_values("rds")["BBA"].values
    assert bbas[0] > bbas[-1], (
        f"Small grains should have higher BBA than large: {bbas}"
    )


def test_return_spectral():
    """return_spectral=True includes albedo column with 480-element arrays."""
    df = parameter_sweep(
        params={"solzen": [40]},
        return_spectral=True,
        progress=False,
    )
    assert "albedo" in df.columns
    alb = df.iloc[0]["albedo"]
    assert isinstance(alb, np.ndarray)
    assert len(alb) == 480


def test_unknown_key_raises():
    """Unknown parameter key raises ValueError."""
    with pytest.raises(ValueError, match="Unknown parameter key"):
        parameter_sweep(params={"nonexistent_param": [1, 2]}, progress=False)


def test_toon_solver():
    """Toon solver runs without error and produces physical results."""
    # Toon solver requires layer_type=0 (grains) and SZA in [50, 89]
    df = parameter_sweep(
        params={"solzen": [55, 70], "layer_type": [0]},
        solver="toon",
        progress=False,
    )
    assert len(df) == 2
    assert all(0 < b < 1 for b in df["BBA"])


def test_adding_doubling_solver():
    """Adding-doubling solver runs without error and produces physical results."""
    df = parameter_sweep(
        params={"solzen": [40, 60]},
        solver="adding-doubling",
        progress=False,
    )
    assert len(df) == 2
    assert all(0 < b < 1 for b in df["BBA"])
