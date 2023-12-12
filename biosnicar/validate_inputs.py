#!/usr/bin/python
import numpy as np


def validate_inputs(ice, illumination, impurities):
    """
    Checks all config to make sure all inputs are valid by calling out to:
    validate_snow_algae()
    validate_ilumination()
    validate_ice()


    Args:
        ice: instance of Ice class
        rt_config: instance of RTConfig class
        model_config: instance of ModelConfig class
        illumination: instance of Illumination class
        impurities: list of instances of Impurity class

    Returns:
        None
    """
    print("\n** Validating model configuration **")
    validate_snow_algae(impurities)
    validate_illumination(illumination)
    validate_ice(ice)
    validate_glacier_algae(impurities)
    return


def validate_snow_algae(impurities):
    """
    Validates snow algae configuration.
    This includes warning the user that OPs are from literature and not yet field
    validated.
    Most importantly, checks that the units of the absorption coefficient and the cell
    concentration match up (m2/cell -> cells/mL, m2/mg -> ppb)

    Args:
        impurities: a list of Impurity objects

    Returns:
        None

    Raises:
        ValueError when units mismatch

    """

    for i, impurity in enumerate(impurities):
        if impurity.file == "snw_alg_r025um_chla020_chlb025_cara150_carb140.nc" and (
            np.sum(impurity.conc) > 0
        ):
            if impurity.unit != 0:
                raise ValueError(
                    "\nyour chosen snow algae file has its absorption cross section in m2/kg, so please express concentration in ppb"
                )
            else:
                print("\n*** WARNING ***")
                print("You have included snow algae as an impurity")
                print("be warned the snow algae optical properties")
                print("are theoretical and yet to be field validated")
        elif (
            impurity.file == "Data/OP_data/480band/lap/SA_Chevrollier2022_r8.99.nc"
            and (np.sum(impurity.conc) > 0)
        ):
            if impurity.unit != 1:
                raise ValueError(
                    "\nyour chosen snow algae file has its absorption cross section in m2/cell, so please express concentration in cells/mL"
                )

    print("snow algae OK")

    return


def validate_glacier_algae(impurities):
    """
    Validates snow algae configuration.
    This includes warning the user that OPs are from literature and not yet field validated.
    Most importantly, checks that the units of the absorption coefficient and the cell
    concentration match up (m2/cell -> cells/mL, m2/mg -> ppb)

    Args:
        impurities: a list of Impurity objects

    Returns:
        None

    Raises:
        ValueError when units mismatch

    """

    for i, impurity in enumerate(impurities):
        if impurity.file == "GA_Chevrollier2022_r4.9_L18.8.nc" and (
            np.sum(impurity.conc) > 0
        ):
            if impurity.unit != 1:
                raise ValueError(
                    "\nyour chosen glacier algae file has its absorption coefficient in m2/cell, so please express concentration in cells/mL"
                )
        elif impurity.file == "Cook2020_glacier_algae_4_40.nc" and (
            np.sum(impurity.conc) > 0
        ):
            if impurity.unit != 0:
                raise ValueError(
                    "\nyour chosen glacier algae file has its absorption coefficient in m2/mg, so please express concentration in ppb"
                )
            else:
                print("\n*** WARNING ***")
                print("You have included glacier algae as an impurity")
                print("be warned the snow algae optical properties")
                print("are theoretical and yet to be field validated")

    print("glacier algae OK")

    return


def validate_illumination(illumination):
    """Validates illumination.

    Args:
        illumination: a list of Impurity objects

    Returns:
        None

    Raises:
        ValueError when SZA or nbr_wvl outside valid range

    """

    if (illumination.solzen > 89) or (illumination.solzen < 1):
        raise ValueError("SZA outside valid range")

    if illumination.nbr_wvl != 480:
        raise ValueError("Illumination has incorrect number of wavelengths")

    if illumination.direct > 1 or illumination.direct < 0:
        raise ValueError("Beam type is incorrect: it should be 0 or 1")

    print("illumination OK")

    return


def validate_ice(ice):
    """Validates ice configuration.

    Args:
        ice: a class containing ice physical constants

    Returns:
        None

    Raises:
        ValueError when lengths of variables are not equal

    """

    equal_len_fields_ice = [
        len(ice.dz),
        len(ice.layer_type),
        len(ice.cdom),
        len(ice.rho),
        len(ice.rds),
        len(ice.shp),
        len(ice.water),
        len(ice.hex_side),
        len(ice.hex_length),
        len(ice.ar),
        ice.nbr_lyr,
    ]

    if not all(length == equal_len_fields_ice[0] for length in equal_len_fields_ice):
        raise ValueError("variables in ice do not have equal lengths")
    for i in range(ice.nbr_lyr):
        if ice.rf != 2 and ice.rds[i] > 1500 and ice.layer_type[i] == 0:
            raise ValueError(
                "Grain size only available up to 1500um with selected ref index"
            )

    print("ice OK")

    return


def validate_model_config(model_config):
    """Validates model configuration.

    Args:
        model_config: a class containing model config variables

    Returns:
        None

    Raises:
        ValueError when wavelengths are incorrect

    """
    if len(model_config.wavelengths) != 480:
        raise ValueError("wavelength range incorrectly configured")
    if model_config.nbr_wvl != 480:
        raise ValueError("nbr_wvl does not equal 480")

    return


if __name__ == "__main__":
    pass
