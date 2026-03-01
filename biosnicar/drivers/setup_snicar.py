#!/usr/bin/python

from pathlib import Path

from biosnicar.utils.load_inputs import load_inputs
from biosnicar.classes import (
    Ice,
    Illumination,
    Impurity,
    ModelConfig,
    PlotConfig,
    RTConfig,
)


def setup_snicar(input_file):
    """Build all model objects from a YAML configuration file.

    Args:
        input_file: Path to the YAML config, or ``"default"`` to use the
            bundled ``inputs.yaml``.

    Returns:
        Tuple of (ice, illumination, rt_config, model_config, plot_config,
        impurities).
    """
    # define input file
    if input_file == "default":
        BIOSNICAR_SRC_PATH = Path(__file__).resolve().parent
        input_file = BIOSNICAR_SRC_PATH.joinpath("../inputs.yaml").as_posix()

    else:
        input_file = input_file
    impurities = build_impurities_array(input_file)
    (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
    ) = build_classes(input_file)

    return (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
        impurities,
    )


def build_classes(input_file):
    """Instantiate model classes from a YAML config file.

    Args:
        input_file: Path to the YAML configuration file.

    Returns:
        Tuple of (ice, illumination, rt_config, model_config, plot_config).
    """

    ice = Ice(input_file)
    illumination = Illumination(input_file)
    rt_config = RTConfig(input_file)
    model_config = ModelConfig(input_file)
    plot_config = PlotConfig(input_file)

    return ice, illumination, rt_config, model_config, plot_config


def build_impurities_array(input_file):
    """Create a list of Impurity instances from YAML config.

    Args:
        input_file: Path to the YAML configuration file.

    Returns:
        List of :class:`~biosnicar.classes.impurity.Impurity` instances.
    """

    inputs = load_inputs(input_file)

    impurities = []

    for i, id in enumerate(inputs["IMPURITIES"]):
        name = inputs["IMPURITIES"][id]["NAME"]
        file = inputs["IMPURITIES"][id]["FILE"]
        coated = inputs["IMPURITIES"][id]["COATED"]
        unit = inputs["IMPURITIES"][id]["UNIT"]
        conc = inputs["IMPURITIES"][id]["CONC"]
        impurities.append(Impurity(file, coated, unit, name, conc))

    return impurities


if __name__ == "__main__":
    pass
