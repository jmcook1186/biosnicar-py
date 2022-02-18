from classes import *


def validate_inputs(ice, rt_config, model_config, illumination, impurities):

    """
    func checks all config to make sure all inputs are valid
    """

    check_toon_instabilities(rt_config, illumination)
    #check_layer_types(ice, rt_config)
    check_algae(impurities)
    check_sza_range(illumination)

    return "ok"


def check_toon_instabilities(rt_config, illumination):
    """
    the toon solver has instabiolities
    """
    if rt_config.toon and illumination.solzen < 50:  # pylint: disable=W0143
        raise ValueError("INVALID SOLZEN: outside valid range for toon solver")

    return


# def check_layer_types(ice, rt_config):
#     """
#     checks solver and ice layer types match up. If granular and AD just
#     warn user faster option is available.
#     If solid ice and toon raise value error as this is invalid.
#     """
#     if np.sum(ice.layer_type) < 1 and rt_config.add_double:  # pylint: disable=W0143
#         # just warn user but let program continue - in some cases
#         # AD method preferable (stable over complete range of SZA)
#         print("*** WARNING ***")
#         print("No solid ice layers - toon solver is faster")
#         print("Toggle toon=True and add_double=False to use it.\n")

#     if np.sum(ice.layer_type) > 0 and rt_config.toon:
#         raise ValueError("There are ice layers - use the adding-doubling solver")

#     return


def check_algae(impurities):
    """
    If snow algae is present warn user that the OPs are from literature
    """
    snw_alg = False
    for i, impurity in enumerate(impurities):
        if ("snw_alg" in impurity.name) and (np.sum(impurity.conc) > 0):
            snw_alg = True

    if snw_alg:
        print("*** WARNING ***")
        print("You have included snow algae as an impurity")
        print("be warned the snow algae optical properties")
        print("are theoretical and yet to be field validated")

    return


def check_sza_range(illumination):
    if (illumination.solzen > 89) | (illumination.solzen < 1):  # pylint: disable=W0143
        raise ValueError(
            "SZA outside of valid range for AD solver - please choose between 1-89 degrees"
        )

    return
