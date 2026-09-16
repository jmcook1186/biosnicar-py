"""Tests for the run_model() entry point."""

import numpy as np
import pytest

from biosnicar.drivers.run_model import run_model
from biosnicar.classes.outputs import Outputs


def test_default_returns_outputs():
    """Default call returns an Outputs object with a valid BBA."""
    outputs = run_model()
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_solzen_changes_bba():
    """Changing solzen produces a different BBA."""
    out1 = run_model(solzen=30)
    out2 = run_model(solzen=70)
    assert out1.BBA != out2.BBA


def test_float_solzen_does_not_raise():
    """A float solar zenith angle must not raise (regression).

    The solar-flux LUT is keyed by zero-padded integer degrees (SZA00..SZA89),
    while the public API documents solzen as a float. A float (or single-digit)
    zenith must round to the integer key rather than leaking "SZA50.0" into the
    lookup and raising KeyError.
    """
    integer = run_model(solzen=50)
    for sz in (50.0, 49.6, 5, 5.0):
        out = run_model(solzen=sz)
        assert 0 < out.BBA < 1
    # 50.0 and 49.6 both round to 50 -> identical result to the integer call
    assert run_model(solzen=50.0).BBA == integer.BBA
    assert run_model(solzen=49.6).BBA == integer.BBA


def test_override_rds_changes_bba():
    """Changing grain radius produces a different BBA."""
    out1 = run_model(rds=200)
    out2 = run_model(rds=5000)
    assert out1.BBA != out2.BBA
    # Larger grains absorb more -> lower BBA
    assert out1.BBA > out2.BBA


def test_outputs_have_spectral_albedo():
    """Outputs always contain 480-element spectral albedo."""
    outputs = run_model()
    assert hasattr(outputs, "albedo")
    assert len(outputs.albedo) == 480


def test_invalid_solver_raises():
    """Invalid solver name raises ValueError."""
    with pytest.raises(ValueError, match="Unknown solver"):
        run_model(solver="invalid-solver")


def test_invalid_override_key_raises():
    """Unknown override key raises ValueError."""
    with pytest.raises(ValueError, match="Unknown override key"):
        run_model(nonexistent_key=42)


def test_get_albedo_still_works():
    """get_albedo.get() still works as a thin wrapper."""
    from biosnicar.drivers.get_albedo import get

    albedo = get("adding-doubling", plot=False, validate=False)
    assert len(albedo) == 480
    assert all(0 <= a <= 1 for a in albedo)


def test_impurity_override():
    """Impurity concentration override affects BBA."""
    out_clean = run_model()
    out_dirty = run_model(black_carbon=50000)
    # More black carbon -> lower albedo
    assert out_dirty.BBA < out_clean.BBA


def test_biosnicar_module_run_model():
    """biosnicar.run_model() works as a top-level convenience."""
    import biosnicar

    outputs = biosnicar.run_model(solzen=50)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_shp_scalar():
    """Scalar shp override broadcasts to all layers."""
    outputs = run_model(shp=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_grain_ar_scalar():
    """Scalar grain_ar override broadcasts to all layers."""
    outputs = run_model(grain_ar=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_hex_geometry():
    """Hexagonal column geometry overrides are accepted."""
    outputs = run_model(hex_side=5000, hex_length=5000)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_water_coating():
    """Water coating override is accepted."""
    outputs = run_model(water=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_cdom():
    """CDOM override is accepted."""
    outputs = run_model(cdom=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_shp_fctr():
    """Shape factor override is accepted."""
    outputs = run_model(shp_fctr=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_override_multiple_new_keys():
    """Multiple new ice overrides can be combined."""
    outputs = run_model(shp=0, grain_ar=0, cdom=0, water=0, shp_fctr=0)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_layer_count_change_resizes_all():
    """Changing layer count via dz resizes all attributes including impurities."""
    outputs = run_model(
        dz=[0.01, 0.02, 0.03],
        rds=[500, 800, 1000],
        rho=[400, 600, 700],
        shp=[0, 0, 0],
        snow_algae=[10000, 0, 0],
    )
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1


def test_layer_count_change_without_impurity_override():
    """Changing layer count works even when no impurity overrides are given."""
    outputs = run_model(
        dz=[0.02, 0.05, 0.05, 0.05, 0.83],
        rds=[800, 900, 1000, 1100, 1200],
        rho=[500, 600, 700, 750, 800],
    )
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1
