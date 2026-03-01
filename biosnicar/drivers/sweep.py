"""Parameter sweep API for running biosnicar across ranges of inputs.

Example usage::

    from biosnicar.drivers.sweep import parameter_sweep

    df = parameter_sweep(
        params={
            "solzen": [30, 40, 50, 60, 70],
            "rds": [100, 200, 500, 1000],
        }
    )
    df.pivot_table(values="BBA", index="solzen", columns="rds").plot()

"""

import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd

from biosnicar.drivers.setup_snicar import setup_snicar
from biosnicar.optical_properties.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.rt_solvers.adding_doubling_solver import adding_doubling_solver
from biosnicar.rt_solvers.toon_rt_solver import toon_solver

# Regex for impurity concentration keys like "impurity.0.conc"
_IMPURITY_CONC_RE = re.compile(r"^impurity\.(\d+)\.conc$")

# Parameter keys that require recalculating irradiance
_ILLUMINATION_KEYS = {"solzen", "direct", "incoming"}

# Parameter keys that apply to the ice object (broadcast to all layers)
_ICE_BROADCAST_KEYS = {"rds", "rho", "dz", "lwc", "layer_type"}

_VALID_KEYS = _ILLUMINATION_KEYS | _ICE_BROADCAST_KEYS

_OUTPUT_SCALARS = [
    "BBA",
    "BBAVIS",
    "BBANIR",
    "abs_slr_tot",
    "abs_vis_tot",
    "abs_nir_tot",
    "abs_slr_btm",
    "total_insolation",
]


def parameter_sweep(
    params,
    solver="adding-doubling",
    input_file="default",
    return_spectral=False,
    progress=True,
):
    """Run biosnicar over the Cartesian product of parameter values.

    Args:
        params: Dict mapping parameter names to lists of values.
            Supported keys: ``solzen``, ``direct``, ``incoming``, ``rds``,
            ``rho``, ``dz``, ``layer_type``, ``impurity.{i}.conc``.
        solver: ``"adding-doubling"`` (default) or ``"toon"``.
        input_file: Path to YAML config or ``"default"``.
        return_spectral: If True, include an ``albedo`` column with 480-element
            numpy arrays.
        progress: If True, display a tqdm progress bar.

    Returns:
        :class:`pandas.DataFrame` with one row per parameter combination.

    Raises:
        ValueError: If an unrecognised parameter key is supplied.
    """
    _validate_keys(params)

    # Resolve "default" to actual path so it can be reused in recalculations
    if input_file == "default":
        input_file = Path(__file__).resolve().parent.joinpath("../inputs.yaml").as_posix()

    # Build baseline objects once
    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar(
        input_file
    )

    # Resolve solver function
    if solver == "adding-doubling":
        solve = _run_adding_doubling
    elif solver == "toon":
        solve = _run_toon
    else:
        raise ValueError(f"Unknown solver {solver!r}; use 'adding-doubling' or 'toon'")

    # Build Cartesian product
    keys = list(params.keys())
    value_lists = [params[k] for k in keys]
    combos = list(itertools.product(*value_lists))

    # Optional progress bar
    iterator = combos
    if progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(combos, desc="parameter_sweep")
        except ImportError:
            pass

    results = []
    for combo in iterator:
        combo_dict = dict(zip(keys, combo))

        # Apply parameter mutations
        _apply_params(combo_dict, ice, illumination, impurities, input_file)

        # Forward model
        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )
        outputs = solve(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)

        # Collect results
        row = dict(combo_dict)
        for attr in _OUTPUT_SCALARS:
            row[attr] = getattr(outputs, attr)
        if return_spectral:
            row["albedo"] = np.array(outputs.albedo)
        results.append(row)

    return pd.DataFrame(results)


def _validate_keys(params):
    """Raise ValueError for any unrecognised parameter key."""
    for key in params:
        if key in _VALID_KEYS:
            continue
        if _IMPURITY_CONC_RE.match(key):
            continue
        raise ValueError(
            f"Unknown parameter key {key!r}. Supported keys: "
            f"{sorted(_VALID_KEYS)} and 'impurity.{{i}}.conc'."
        )


def _apply_params(combo_dict, ice, illumination, impurities, input_file):
    """Mutate model objects in-place according to a single parameter combination."""
    needs_irradiance = False
    needs_refractive = False

    for key, value in combo_dict.items():
        # Illumination scalars
        if key == "solzen":
            illumination.solzen = value
            needs_irradiance = True
        elif key == "direct":
            illumination.direct = value
            needs_irradiance = True
        elif key == "incoming":
            illumination.incoming = value
            needs_irradiance = True

        # Ice broadcast keys
        elif key in _ICE_BROADCAST_KEYS:
            broadcast = [value] * ice.nbr_lyr
            setattr(ice, key, broadcast)
            needs_refractive = True

        # Impurity concentration
        else:
            m = _IMPURITY_CONC_RE.match(key)
            if m:
                idx = int(m.group(1))
                impurities[idx].conc = [value] * ice.nbr_lyr

    if needs_refractive:
        ice.calculate_refractive_index(input_file)
    if needs_irradiance:
        illumination.calculate_irradiance()


def _run_adding_doubling(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config):
    return adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config)


def _run_toon(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config):
    return toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config)
