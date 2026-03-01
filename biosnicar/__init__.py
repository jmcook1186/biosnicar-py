from pathlib import Path

__version__ = "2.1.0"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"


def run_model(*args, **kwargs):
    """Run the BioSNICAR forward model. See :func:`biosnicar.drivers.run_model.run_model`."""
    from biosnicar.drivers.run_model import run_model as _run_model

    return _run_model(*args, **kwargs)
