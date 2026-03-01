"""Single entry point for running the BioSNICAR forward model.

Provides :func:`run_model`, which sets up the model from a YAML config,
applies optional keyword overrides, runs the radiative transfer solver,
and returns an :class:`~biosnicar.classes.outputs.Outputs` object.

Example::

    from biosnicar.drivers.run_model import run_model

    outputs = run_model(solzen=50, rds=1000)
    print(outputs.BBA)
"""

import re

from biosnicar.drivers.setup_snicar import setup_snicar
from biosnicar.optical_properties.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.rt_solvers.adding_doubling_solver import adding_doubling_solver
from biosnicar.rt_solvers.toon_rt_solver import toon_solver
from biosnicar.utils.display import display_out_data, plot_albedo
from biosnicar.utils.validate_inputs import validate_inputs

# Regex for impurity concentration keys like "impurity_0_conc"
_IMPURITY_CONC_RE = re.compile(r"^impurity_(\d+)_conc$")

# Parameter keys that require recalculating irradiance
_ILLUMINATION_KEYS = {"solzen", "direct", "incoming"}

# Parameter keys that apply to the ice object (broadcast to all layers)
_ICE_BROADCAST_KEYS = {"rds", "rho", "dz", "lwc", "layer_type"}


def run_model(
    input_file="default",
    solver="adding-doubling",
    validate=False,
    plot=False,
    **overrides,
):
    """Run the BioSNICAR forward model and return outputs.

    Builds model objects from *input_file*, applies any keyword *overrides*
    to ice/illumination/impurity parameters, then runs the full pipeline
    (optical properties -> impurity mixing -> radiative transfer).

    Args:
        input_file: Path to YAML config or ``"default"`` for the bundled
            ``inputs.yaml``.
        solver: ``"adding-doubling"`` (default) or ``"toon"``.
        validate: If True, run input validation before the forward model.
        plot: If True, display a spectral albedo plot after the run.
        **overrides: Parameter overrides applied before running the model.
            Supported keys:

            - **solzen** (*float*) — solar zenith angle (degrees)
            - **direct** (*int*) — 1 for direct beam, 0 for diffuse
            - **incoming** (*int*) — irradiance spectrum index (0–6)
            - **rds** (*int | list*) — grain/bubble radius (µm), broadcast
              to all layers if scalar
            - **rho** (*int | list*) — layer density (kg/m³)
            - **dz** (*float | list*) — layer thickness (m)
            - **lwc** (*float | list*) — liquid water content
            - **layer_type** (*int | list*) — 0=grains, 1=solid ice, etc.
            - **impurity_{i}_conc** (*float | list*) — concentration for
              impurity *i* (0-based index), broadcast to ``[value, 0, ...]``
              if scalar

    Returns:
        :class:`~biosnicar.classes.outputs.Outputs` with albedo, BBA,
        absorbed fluxes, etc.

    Raises:
        ValueError: If *solver* is not recognised or an override key is
            unknown.
    """
    # Resolve "default" to actual path so it can be reused in recalculations
    if input_file == "default":
        from pathlib import Path

        input_file = Path(__file__).resolve().parent.joinpath(
            "../inputs.yaml"
        ).as_posix()

    ice, illumination, rt_config, model_config, plot_config, impurities = (
        setup_snicar(input_file)
    )

    if overrides:
        _apply_overrides(overrides, ice, illumination, impurities, input_file)

    if validate:
        validate_inputs(ice, illumination, impurities)

    # Optical properties
    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(
        ssa_snw, g_snw, mac_snw, ice, impurities, model_config
    )

    # Radiative transfer
    if solver == "adding-doubling":
        outputs = adding_doubling_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config
        )
    elif solver == "toon":
        outputs = toon_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )
    else:
        raise ValueError(
            f"Unknown solver {solver!r}; use 'adding-doubling' or 'toon'"
        )

    if plot:
        plot_albedo(plot_config, model_config, outputs.albedo)
    display_out_data(outputs)

    return outputs


def _apply_overrides(overrides, ice, illumination, impurities, input_file):
    """Mutate model objects in-place according to keyword overrides."""
    needs_irradiance = False
    needs_refractive = False

    for key, value in overrides.items():
        # Illumination scalars
        if key in _ILLUMINATION_KEYS:
            setattr(illumination, key, value)
            needs_irradiance = True

        # Ice broadcast keys
        elif key in _ICE_BROADCAST_KEYS:
            if isinstance(value, list):
                setattr(ice, key, value)
                ice.nbr_lyr = len(ice.dz)
            else:
                setattr(ice, key, [value] * ice.nbr_lyr)
            needs_refractive = True

        # Impurity concentration (impurity_0_conc, impurity_1_conc, ...)
        else:
            m = _IMPURITY_CONC_RE.match(key)
            if m:
                idx = int(m.group(1))
                if isinstance(value, list):
                    impurities[idx].conc = value
                else:
                    # Broadcast: put concentration in first layer, zero elsewhere
                    impurities[idx].conc = [value] + [0] * (ice.nbr_lyr - 1)
            else:
                raise ValueError(
                    f"Unknown override key {key!r}. Supported: "
                    f"{sorted(_ILLUMINATION_KEYS | _ICE_BROADCAST_KEYS)} "
                    f"and 'impurity_{{i}}_conc'."
                )

    if needs_refractive:
        ice.calculate_refractive_index(input_file)
    if needs_irradiance:
        illumination.calculate_irradiance()
