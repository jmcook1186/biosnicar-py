#################################################################
# this file will be deleted - using as PoC for new testing syntax
#################################################################

from setup_snicar import *
from classes import *
from column_OPs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver
from validate_inputs import *
from plot import *


# first setup snicar with vals in config files ("defaults")
ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
status = validate_inputs(ice, rt_config, model_config, illumination, impurities)

# array of vars to iterate over
densList = [400, 500, 600, 700, 800]
reffList = [200, 400, 600, 800, 1000]

# set up loop
for density in densList:
    for reff in reffList:

        # update fields in existign class
        ice.rho = [density] * len(densList)
        ice.rds = [reff] * len(reffList)

        # OP calcs
        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )

        # get albedo
        outputs = adding_doubling_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

        print(outputs.BBA)
