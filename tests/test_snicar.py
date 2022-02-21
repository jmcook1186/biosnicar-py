import sys

from src.classes import Impurity

# make sure we can import from/src
sys.path.append("./src")
from setup_snicar import *
from classes import *
from column_OPs import *
from toon_rt_solver import toon_solver
from adding_doubling_solver import adding_doubling_solver
import random
import matplotlib.pyplot as plt
import pytest


"""
To run configure these tests, update the values in conftest.py
Then navigate to the tests folder and run

`pytest`

The tests will automatically run - green dots indicate tests
passing successfully. A plot of N random spectra pairs will
be saved to the /tests folder.

The fuzzer exists to run snocar with a wide range of input variables
to check that no combinations break the code. It is quite memory intensive
to fuzz over a very large parameter space. 10^3 runs is ok on a decent
spec laptop. I have divided the fuzzer into two separate functions. One
has coverage for "config" variables that set up the radiative transfer
e.g. direct vs diffuse, approximation type, etc. The other is more for
conditions of the ice column, e.g. density, effective radius, LAPs.

To toggle the fuzzer on/off change the value of "fuzz" in conftest.py

"""

def test_v3(new_benchmark_toon):
    if new_benchmark_toon:
        print("generating benchmark data using params equivalent to snicarv3 (Toon solver)")
        # first setup snicar with vals in config files
        ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
        status = validate_inputs(ice, rt_config, model_config, illumination, impurities)
        ice, impurities, illumination, rt_config, model_config = match_matlab_config(ice, illumination, rt_config, model_config)


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
                                tau, ssa, g, L_snw, ice, illumination, model_config
                                )

                                specOut[counter, 0:480] = outputs.albedo
                                specOut[counter, 480] = outputs.BBA
                                counter +=1


        np.savetxt("./tests/test_data/py_benchmark_data_toon.csv", specOut, delimiter=",")

    else:
        pass

    return


def test_v4(new_benchmark_ad):
    
    if new_benchmark_ad:
        ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
        ice, illumination, impurities, rt_config, model_config = match_matlab_config(ice, illumination, rt_config, model_config)
        
        print("generating benchmark data using params equivalent to snicarv4 (AD solver)")

        lyrList = [0,1]
        densList = [400, 500, 600, 700, 800]
        reffList = [200, 400, 600, 800, 1000]
        zenList = [30, 40, 50, 60]
        bcList = [500, 1000, 2000]
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
        )

        assert(ncols==3000)

        specOut = np.zeros(shape=(ncols, 481))
        counter = 0
        for layer_type in lyrList:
            for density in densList:
                for reff in reffList:
                    for zen in zenList:
                        for bc in bcList:
                            for dz in dzList:
                                
                                ice.dz = dz
                                ice.nbr_lyr = 5
                                ice.layer_type = [layer_type]*len(ice.dz)
                                ice.rho = [density]*len(ice.dz)
                                ice.rds = [reff]*len(ice.dz)
                                illumination.solzen = zen
                                impurities[0].conc = [bc]*len(ice.dz) #bc in all layers
                                
                                assert(len(impurities[0].conc)==5)
                                assert(impurities[0].name=='bc')
                                assert(impurities[0].cfactor==1)
                                assert(ice.nbr_lyr==5)
                                assert(len(impurities)==1)
                                assert(impurities[0].conc==[bc, bc, bc, bc, bc])
                                assert(illumination.incoming==4)
                                assert(illumination.solzen == zen)
                                assert(ice.shp == [0, 0, 0, 0, 0])
                                assert(ice.shp_fctr == [0, 0, 0, 0, 0])

                                ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
                                tau, ssa, g, L_snw = mix_in_impurities(
                                    ssa_snw, g_snw, mac_snw, ice, impurities, model_config
                                )

                                outputs = adding_doubling_solver(
                                tau, ssa, g, L_snw, ice, illumination, model_config
                                )

                                specOut[counter, 0:480] = outputs.albedo
                                specOut[counter, 480] = outputs.BBA
                                counter +=1

        np.savetxt("./tests/test_data/py_benchmark_data.csv", specOut, delimiter=",")
    
    else:
        pass

    return



def test_realistic_BBA(get_matlab_data, get_python_data):
    # are the values predicted by the model always physical (i.e. between 0-1)
    # do the files have the right shape and size?

    mat = get_matlab_data
    py = get_python_data

    bb_py = py.loc[:, 481]
    bb_mat = mat.loc[:, 481]

    assert len(bb_py) == len(bb_mat)
    assert bb_py[bb_py > 1].count() == 0 and bb_py[bb_py < 0].count() == 0
    assert bb_mat[bb_mat > 1].count() == 0 and bb_mat[bb_mat < 0].count() == 0


def test_compare_pyBBA_to_matBBA(get_matlab_data, get_python_data, set_tolerance):
    # check the BBA predicted for each run matches to within tolerance between
    # the two models
    mat = get_matlab_data
    py = get_python_data
    tol = set_tolerance
    bb_py = py.loc[:, 481]
    bb_mat = mat.loc[:, 481]
    error = np.array(abs(bb_mat - bb_py))
    assert len(error[error > tol]) == 0


def match_matlab_config(ice, illumination, rt_config, model_config):
    
    nbr_lyr = 5
    # make sure ice config matches matlab benchmark
    ice.ri = 2
    ice.shp = [0]*nbr_lyr
    ice.shp_fctr = [0]*nbr_lyr
    ice.ar = [0]*nbr_lyr
    ice.sfc = np.array([0.25]*model_config.nbr_wvl)
    ice.cdom = [0]*nbr_lyr
    ice.water = [0]*nbr_lyr
    ice.nbr_lyr = nbr_lyr
    ice.layer_type = [0]*nbr_lyr

    illumination.incoming = 4
    illumination.direct = 1

    impurities = []

    # make sure impurities[0] is bc
    # (same bc used by matlab model)

    conc = [0]*nbr_lyr
    impurity0 = Impurity(model_config.dir_base, "mie_sot_ChC90_dns_1317.nc", False, 1, 0, "bc", conc)
    impurities.append(impurity0)

    assert((impurities[0].name=="bc") and (impurities[0].file == "mie_sot_ChC90_dns_1317.nc"))

    return ice, illumination, impurities, rt_config, model_config


def test_compare_pySPEC_to_matSPEC(get_matlab_data, get_python_data, set_tolerance):
    # check that for each individual wavelenght, the spectral albedo
    # matches to within tolerance between the two models
    mat = get_matlab_data
    py = get_python_data
    tol = set_tolerance
    bb_spec = py.loc[:, :480]
    bb_spec = mat.loc[:, :480]
    error = np.array(abs(bb_spec - bb_spec))
    assert len(error[error > tol]) == 0


def test_compare_pySPEC_to_matSPEC_TOON(
    get_matlab_data_toon, get_python_data_toon, set_tolerance
):
    # check that for each individual wavelenght, the spectral albedo
    # matches to within tolerance between the two models
    mat = get_matlab_data_toon
    py = get_python_data_toon
    tol = set_tolerance
    bb_spec = py.loc[:, :480]
    bb_spec = mat.loc[:, :480]
    error = np.array(abs(bb_spec - bb_spec))
    assert len(error[error > tol]) == 0


def test_plot_random_spectra_pairs(get_matlab_data, get_python_data, get_n_spectra):
    # grabs n spectra and plots the python and matlab versions
    # for visual comparison
    mat = get_matlab_data
    py = get_python_data
    n_spec = get_n_spectra
    idxs = random.sample(range(0, py.shape[0]), n_spec)

    wavelength = np.arange(200, 5000, 10)
    py_spectra = py.iloc[0:-1, idxs]
    mat_spectra = mat.iloc[0:-1, idxs]

    plt.plot(wavelength, py_spectra)
    plt.plot(wavelength, mat_spectra, linestyle=None, marker="x")
    plt.xlabel("wavelength (nm)")
    plt.ylabel("Albedo")
    plt.title("solid lines = Python\ncrosses = Matlab")

    plt.savefig("./tests/test_data/py_mat_comparison.png")


# each parametrize decorator defines a range of values for
# a specific input parameter
@pytest.mark.parametrize("direct", [0, 1])
@pytest.mark.parametrize("aprx_typ", [1, 2, 3])
@pytest.mark.parametrize("incoming", [0, 1, 2, 3, 4, 5, 6])
@pytest.mark.parametrize("rf", [0, 1, 2])
def test_config_fuzzer(direct, aprx_typ, incoming, rf, fuzz):
    """
    ensures code runs and gives valid BBA with combinations of input vals
    runs for both solvers
    """
    
    if fuzz:
        ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
        status = validate_inputs(ice, rt_config, model_config, illumination, impurities)
        ice, impurities, illumination, rt_config, model_config = match_matlab_config(ice, illumination, rt_config, model_config)


        rt_config.aprx_typ = aprx_typ
        illumination.direct = direct
        illumination.incoming = incoming
        ice.rf = rf 

        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )
        outputs_toon = toon_solver(
        tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

        outputs_ad = toon_solver(
        tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

    else:
        pass

    return


@pytest.mark.parametrize("rds", [1000, 5000, 10000])
@pytest.mark.parametrize("rho", [400, 600, 800])
@pytest.mark.parametrize("solzen", [30, 50, 70])
@pytest.mark.parametrize("cfactor", [1, 30])
@pytest.mark.parametrize("dust", [0, 50000])
@pytest.mark.parametrize("algae", [0, 50000])
def test_var_fuzzer(rds, rho, solzen, cfactor, dust, algae, fuzz):
    """
    ensures code runs and gives valid BBA with combinations of input vals
    """

    if fuzz:
        ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar()
        status = validate_inputs(ice, rt_config, model_config, illumination, impurities)
        ice, impurities, illumination, rt_config, model_config = match_matlab_config(ice, illumination, rt_config, model_config)

        impurities = []

        conc1 = [0]*len(ice.dz)
        conc1[0] = algae
        impurity0 = Impurity(model_config.dir_base, "mie_sot_ChC90_dns_1317.nc", False, cfactor, 0, "bc", conc1)
        impurities.append(impurity0)

        conc2 = [0]*len(ice.dz)
        conc2[0] = dust
        impurity1 = Impurity(model_config.dir_base, "dust_balkanski_central_size1.nc", False, cfactor, 0, "dust1", conc2)
        impurities.append(impurity1)


        ice.rds = [rds]*len(ice.dz)
        ice.rho = [rho]*len(ice.dz)
        illumination.solzen = solzen


        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, impurities, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )
        outputs_toon = toon_solver(
        tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

        outputs_ad = toon_solver(
        tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

    else:
        pass


    return



