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
        name = inputs["IMPURITIES"][id]["name"]
        file = inputs["IMPURITIES"][id]["file"]
        cfactor = inputs["IMPURITIES"][id]["cfactor"]
        coated = inputs["IMPURITIES"][id]["coated"]
        unit = inputs["IMPURITIES"][id]["unit"]
        conc = inputs["IMPURITIES"][id]["conc"]
        impurities.append(Impurity(dir_base, file, coated, cfactor, unit, name, conc))

    return impurities

