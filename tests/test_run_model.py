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
    out_dirty = run_model(impurity_0_conc=50000)
    # More black carbon -> lower albedo
    assert out_dirty.BBA < out_clean.BBA


def test_biosnicar_module_run_model():
    """biosnicar.run_model() works as a top-level convenience."""
    import biosnicar

    outputs = biosnicar.run_model(solzen=50)
    assert isinstance(outputs, Outputs)
    assert 0 < outputs.BBA < 1
