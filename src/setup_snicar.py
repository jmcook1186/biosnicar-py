import yaml
from classes import *
import csv


def setup_snicar():

    impurities = build_impurities_array()
    ice, illumination, rt_config, model_config, plot_config = build_classes()

    return ice, illumination, rt_config, model_config, plot_config, impurities


def build_classes():

    ice = Ice()
    illumination = Illumination()
    rt_config = RTConfig()
    model_config = ModelConfig()
    plot_config = PlotConfig()

    return ice, illumination, rt_config, model_config, plot_config


def build_impurities_array():

    """
    creates an array of impurities - each one an instance of Impurity with
    properties defined in impurity_config.yaml
    """

    with open("./src/inputs.yaml", "r") as ymlfile:
        inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

    impurities = []
    dir_base = inputs["PATHS"]["DIR_BASE"]

    for i, id in enumerate(inputs["IMPURITIES"]):
        name = inputs["IMPURITIES"][id]["NAME"]
        file = inputs["IMPURITIES"][id]["FILE"]
        cfactor = inputs["IMPURITIES"][id]["CFACTOR"]
        coated = inputs["IMPURITIES"][id]["COATED"]
        unit = inputs["IMPURITIES"][id]["UNIT"]
        conc = inputs["IMPURITIES"][id]["CONC"]
        impurities.append(Impurity(dir_base, file, coated, cfactor, unit, name, conc))

    return impurities
