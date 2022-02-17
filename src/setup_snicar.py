import yaml
from classes import *
import csv


def get_config():
    """

    get_config() creates config objects from the relevant .yaml files to

    """
    with open(
        "/home/joe/Code/BioSNICAR_GO_PY/src/impurity_config.yaml", "r"
    ) as ymlfile:
        impurity_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    with open("/home/joe/Code/BioSNICAR_GO_PY/src/rtm_config.yaml", "r") as ymlfile:
        rtm_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    with open(
        "/home/joe/Code/BioSNICAR_GO_PY/src/ice_physical_config.yaml", "r"
    ) as ymlfile:
        ice_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    with open("/home/joe/Code/BioSNICAR_GO_PY/src/model_config.yaml", "r") as ymlfile:
        model_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    return impurity_cfg, rtm_cfg, ice_cfg, model_cfg


def get_wavelengths(model_cfg):
    """
    reads in wavelengths from text file
    """
    with open(model_cfg["PATHS"]["WVL"], newline="") as f:
        reader = csv.reader(f)
        wvl = list(reader)

    return wvl


def build_impurities_array():
    """
    creates an array of impurities - each one an instance of Impurity with
    properties defined in impurity_config.yaml
    """

    with open("/home/joe/Code/BioSNICAR_GO_PY/src/impurity_config.yaml", "r") as ymlfile:
        impurity_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
    
    with open("/home/joe/Code/BioSNICAR_GO_PY/src/model_config.yaml", "r") as ymlfile:
        model_cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    impurities = []
    dir_base = model_cfg["PATHS"]["DIR_BASE"]
    for i, id in enumerate(impurity_cfg["impurities"]):
        name = impurity_cfg["impurities"][id]["name"]
        file = impurity_cfg["impurities"][id]["file"]
        cfactor = impurity_cfg["impurities"][id]["cfactor"]
        coated = impurity_cfg["impurities"][id]["coated"]
        unit = impurity_cfg["impurities"][id]["unit"]
        conc = impurity_cfg["impurities"][id]["conc"]
        impurities.append(Impurity(dir_base, file, coated, cfactor, unit, name, conc))

    return impurities


def build_ice_column(rtm_cfg, model_cfg):
    """
    creates an instance of Ice with properties from ice_physical_config.yaml
    """
    return Ice(rtm_cfg, model_cfg)


def get_illumination(rtm_cfg):
    illumination = Illumination(rtm_cfg)
    illumination.calculate_flx_slr()
    return illumination