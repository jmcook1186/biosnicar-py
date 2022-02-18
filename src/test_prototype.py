#################################################################
# this file will be deleted - using as PoC for new testing syntax
#################################################################

from setup_snicar import *
from classes import *
from column_OPs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver
from validate_inputs import *

def test_v3():
    
    print("generating benchmark data using params equivalent to snicarv3 (Toon solver)")
    # first setup snicar with vals in config files ("defaults")
    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
    status = validate_inputs(ice, rt_config, model_config, illumination, impurities)
    ice, impurities, illumination, rt_config, model_config = match_matlab_config(ice, impurities, illumination, rt_config, model_config)

    lyrList = [0]
    densList = [400, 500, 600, 700, 800]
    reffList = [200, 400, 600, 800, 1000]
    zenList = [30, 40, 50, 60, 70]
    bcList = [500, 1000, 1500, 2000]
    dust1 = [0]
    dust5 = [0]
    dzList = [
        [0.02, 0.04, 0.06, 0.08, 0.1],
        [0.04, 0.06, 0.08, 0.10, 0.15],
        [0.05, 0.10, 0.15, 0.2, 0.5],
        [0.15, 0.2, 0.25, 0.3, 0.5],
        [0.5, 0.5, 0.5, 1, 10],
    ]
    
    ncols = (
    len(lyrList)
    * len(densList)
    * len(reffList)
    * len(zenList)
    * len(bcList)
    * len(dzList)
    * len(dust1)
    * len(dust5)
    )

    specOut = np.zeros(shape=(ncols, 481))
    counter = 0

    for layer_type in lyrList:
        for density in densList:
            for reff in reffList:
                for zen in zenList:
                    for bc in bcList:
                        for dz in dzList:

                            ice.layer_type = [layer_type]*len(ice.layer_type)
                            ice.rho = [density]*len(ice.rho)
                            ice.rds = [reff]*len(ice.rds)
                            illumination.sza = zen
                            impurities[0].conc[0] = bc
                            ice.dz = dz
                            
                            ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
                            tau, ssa, g, L_snw = mix_in_impurities(
                                ssa_snw, g_snw, mac_snw, ice, impurities, model_config
                            )
                            outputs = toon_solver(
                            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
                            )

                            specOut[counter, 0:480] = outputs.albedo
                            specOut[counter, 480] = outputs.BBA
                            counter +=1


    np.savetxt("py_benchmark_data.csv", specOut, delimiter=",")

    return


def test_v4():
    
    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar(defaults)
    status = validate_inputs(ice, rt_config, model_config, illumination, impurities)
    ice, impurities, illumination, rt_config, model_config = match_matlab_config(ice, impurities, illumination, rt_config, model_config)
    
    print("generating benchmark data using params equivalent to snicarv4 (AD solver)")

    lyrList = [0, 1]
    densList = [400, 500, 600, 700, 800]
    reffList = [200, 400, 600, 800, 1000]
    zenList = [30, 40, 50, 60]
    bcList = [0, 1000, 2000]
    dust1List = [0, 10000, 20000, 50000]
    dust5List = [0, 10000, 20000, 50000]
    dzList = [
        [0.02, 0.04, 0.06, 0.08, 0.1],
        [0.04, 0.06, 0.08, 0.10, 0.15],
        [0.05, 0.10, 0.15, 0.2, 0.5],
        [0.15, 0.2, 0.25, 0.3, 0.5],
        [0.5, 0.5, 0.5, 1, 10],
    ]

    specOut = np.zeros(shape=(ncols, 481))
    counter = 0
    for layer_type in lyrList:
        for density in densList:
            for reff in reffList:
                for zen in zenList:
                    for bc in bcList:
                        for dz in dzList:

                            ice.layer_type = [layer_type]*len(dz)
                            ice.rho = [density]*len(dz)
                            ice.rds = [reff]*len(ice.rds)
                            illumination.sza = zen
                            impurities[0].conc[0] = bc
                            ice.dz = dz

                            outputs = toon_solver(
                            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
                            )

                            specOut[counter, 0:480] = outputs.albedo
                            specOut[counter, 480] = outputs.BBA
                            counter +=1


    np.savetxt("py_benchmark_data.csv", specOut, delimiter=",")

    return



def match_matlab_config(ice, impurities, illumination, rt_config, model_config):

    nbr_lyr = 5
    # make sure ice config matches matlab benchmark
    ice.ri = 2
    ice.shp = [0]*nbr_lyr
    ice.shp_fctr = [0]*nbr_lyr
    ice.ar = [0]*nbr_lyr
    ice.sfc = np.array([0.25]*model_config.nbr_wvl)

    illumination.incoming = 4
    illumination.direct = 1

    # make sure impurities[0] is bc
    # (same bc used by matlab model)
    conc = [0]
    impurity = Impurity(model_config.dir_base, "mie_sot_ChC90_dns_1317.nc", False, 1, 0, "bc", conc)
    impurities[0] = impurity


    return ice, impurities, illumination, rt_config, model_config


test_v3()