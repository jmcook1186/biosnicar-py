#!/usr/bin/python
import sys
sys.path.append("./src")
import yaml
from classes import *



def setup_snicar():
    """Builds impurity array and instances of all classes according to config in yaml file.

    Args:
        None

    Returns:
        ice: instance of Ice class
        illumination: instance of Illumination class
        rt_config: instance of RTConfig class
        model_config: instance of ModelConfig class
        plot_config: instance of PlotConfig class
        display_config: instance of DisplayConfig class

    """

    impurities = build_impurities_array()
    (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
    ) = build_classes()

    return (
        ice,
        illumination,
        rt_config,
        model_config,
        plot_config,
        impurities,
    )


def build_classes():
    """Instantiates classes according to config in yaml file.

    Args:
        None

    Returns:
        ice: instance of Ice class
        illumination: instance of Illumination class
        rt_config: instance of RTConfig class
        model_config: instance of ModelConfig class
        plot_config: instance of PlotConfig class
        display_config: instance of DisplayConfig class
    """

    ice = Ice()
    illumination = Illumination()
    rt_config = RTConfig()
    model_config = ModelConfig()
    plot_config = PlotConfig()

    return ice, illumination, rt_config, model_config, plot_config


def build_impurities_array():

    """Creates an array of instances of Impurity.

    creates an array of impurities - each one an instance of Impurity with
    properties defined in yaml file.

    Args:
        None

    Returns:
        impurities: array of instances of Impurity
    """

    with open("./src/inputs.yaml", "r") as ymlfile:
        inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

    impurities = []

    for i, id in enumerate(inputs["IMPURITIES"]):
        name = inputs["IMPURITIES"][id]["NAME"]
        file = inputs["IMPURITIES"][id]["FILE"]
        cfactor = inputs["IMPURITIES"][id]["CFACTOR"]
        coated = inputs["IMPURITIES"][id]["COATED"]
        unit = inputs["IMPURITIES"][id]["UNIT"]
        conc = inputs["IMPURITIES"][id]["CONC"]
        impurities.append(Impurity(file, coated, cfactor, unit, name, conc))

    return impurities

if __name__ == '__main__':
    pass