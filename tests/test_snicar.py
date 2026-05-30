#!/usr/bin/python
"""Runs benchmarking and fuzzing tests on BioSNICAR.

To run configure these tests, update the values in conftest.py
Then navigate to the tests folder and run

`pytest .`

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

import random

import matplotlib.pyplot as plt
import numpy as np
import pytest
from biosnicar.rt_solvers.adding_doubling_solver import adding_doubling_solver
from biosnicar.classes import Impurity
from biosnicar.optical_properties.column_OPs import get_layer_OPs, mix_in_impurities
from biosnicar.drivers.setup_snicar import setup_snicar
from biosnicar.rt_solvers.toon_rt_solver import toon_solver


def test_AD_solver(new_benchmark_ad, input_file):
    """Tests Toon solver against SNICAR_ADv4 benchmark.

    This func generates a new file - py_benchmark_data.csv - that contains
    spectral and broadband albedo simulated by BioSNICAR for a range of input
    configurations. The same set of simulations was also run using a previously
    published version of the SNICAR code written in Matlab by Chloe Whicker at
    University of Michigan and run on the UMich server. This function
    only creates the equivalent dataset using BioSNICAR, it doesn't compare the two.

    Equivalence between the Python and Matlab model configuration is controlled by
    a call to match_matlab_config(). This function can be toggled off by setting
    new_benchmark_ad to False in conftest.py.

    Args:
        new_benchmark_ad: Boolean toggling this function on/off

    Returns:
        None but saves py_benchmark_data.csv to ./tests/test_data/

    """
    if new_benchmark_ad:
        (
            ice,
            illumination,
            rt_config,
            model_config,
            plot_config,
            impurities,
        ) = setup_snicar("default")
        
        # make sure the BH83 data is used as per Matlab's version
        model_config.sphere_ice_path = "data/OP_data/480band/ice_spherical_grains_BH83/"
        model_config.bubbly_ice_path = "data/OP_data/480band/bubbly_ice_files_BH83/"

        ice, illumination, impurities, rt_config, model_config = match_matlab_config(
            ice, illumination, rt_config, model_config, input_file
        )

        lyrList = [0, 1]
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

        assert ncols == 3000

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
                                ice.layer_type = [layer_type] * len(ice.dz)
                                ice.rho = [density] * len(ice.dz)
                                ice.lwc = [0] * len(ice.dz)
                                ice.lwc_pct_bbl = [0] * len(ice.dz)
                                ice.rds = [reff] * len(ice.dz)
                                illumination.solzen = zen
                                illumination.calculate_irradiance()
                                impurities[0].conc = [
                                    bc,
                                    bc,
                                    bc,
                                    bc,
                                    bc,
                                ]  # bc in all layers
                                ice.calculate_refractive_index(input_file)
                                illumination.calculate_irradiance()

                                ssa_snw, g_snw, mac_snw = get_layer_OPs(
                                    ice, model_config
                                )
                                tau, ssa, g, L_snw = mix_in_impurities(
                                    ssa_snw,
                                    g_snw,
                                    mac_snw,
                                    ice,
                                    impurities,
                                    model_config,
                                )
                                outputs = adding_doubling_solver(
                                    tau, ssa, g, L_snw, ice, illumination, model_config
                                )

                                specOut[counter, 0:480] = outputs.albedo
                                specOut[counter, 480] = outputs.BBA
                                counter += 1

        np.savetxt("./tests/test_data/py_benchmark_data.csv", specOut, delimiter=",")

    else:
        pass

    return


def test_AD_solver_clean(new_benchmark_ad_clean, input_file):
    """Tests Toon solver against SNICAR_ADv4 benchmark for impurity-free ice.

    This func generates a new file - py_benchmark_data_clean.csv - that contains
    spectral and broadband albedo simulated by BioSNICAR for a range of input
    configurations. The same set of simulations was also run using a previously
    published version of the SNICAR code written in Matlab by Chloe Whicker at
    University of Michigan and run on the UMich server. This function
    only creates the equivalent dataset using BioSNICAR, it doesn't compare the two.
    The difference between this function and test_v4 is that no impurities are included
    in the model configuration.

    Equivalence between the Python and Matlab model configuration is controlled by
    a call to match_matlab_config(). This function can be toggled off by setting
    new_benchmark_clean to False in conftest.py.

    Args:
        new_benchmark_clean: Boolean toggling this function on/off

    Returns:
        None but saves py_benchmark_data_clean.csv to ./tests/test_data/

    """

    if new_benchmark_ad_clean:
        (
            ice,
            illumination,
            rt_config,
            model_config,
            plot_config,
            impurities,
        ) = setup_snicar("default")
        
        # make sure the BH83 data is used as per Matlab's version
        model_config.sphere_ice_path = "data/OP_data/480band/ice_spherical_grains_BH83/"
        model_config.bubbly_ice_path = "data/OP_data/480band/bubbly_ice_files_BH83/"

        ice, illumination, impurities, rt_config, model_config = match_matlab_config(
            ice, illumination, rt_config, model_config, input_file
        )

        print(
            "generating benchmark data using params equivalent to snicarv4 (AD solver)"
        )

        lyrList = [0, 1]
        densList = [400, 500, 600, 700, 800]
        reffList = [200, 400, 600, 800, 1000]
        zenList = [30, 40, 50, 60]
        bcList = [0]
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
                                ice.layer_type = [layer_type] * len(ice.dz)
                                ice.rho = [density] * len(ice.dz)
                                ice.lwc = [0] * len(ice.dz)
                                ice.lwc_pct_bbl = [0] * len(ice.dz)
                                ice.rds = [reff] * len(ice.dz)
                                illumination.solzen = zen
                                illumination.calculate_irradiance()
                                impurities[0].conc = [
                                    bc,
                                    bc,
                                    bc,
                                    bc,
                                    bc,
                                ]  # bc in all layers
                                ice.calculate_refractive_index(input_file)

                                ssa_snw, g_snw, mac_snw = get_layer_OPs(
                                    ice, model_config
                                )
                                tau, ssa, g, L_snw = mix_in_impurities(
                                    ssa_snw,
                                    g_snw,
                                    mac_snw,
                                    ice,
                                    impurities,
                                    model_config,
                                )

                                outputs = adding_doubling_solver(
                                    tau, ssa, g, L_snw, ice, illumination, model_config
                                )

                                specOut[counter, 0:480] = outputs.albedo
                                specOut[counter, 480] = outputs.BBA
                                counter += 1

        np.savetxt(
            "./tests/test_data/py_benchmark_data_clean.csv", specOut, delimiter=","
        )

    else:
        pass

    return


def test_realistic_bba_ad(get_matlab_data, get_python_data):
    """Tests that BBA values are never >1 or <0.

    Simple test that ensures broadband albedo predicted by model never
    goes outside valid range of 0-1.

    Args:
        get_matlab_data: matlab-snicar-generated csv file of spectral and broadband albedo
        get_python_data: BioSNICAR generated csv file of spectral and broadband albedo

    Returns:
        None

    Raises:
        tests fail if the two datasets differ in length
        tests fail if any BBA values in matlab data are <0 or >1
        tests fail if any BBA values in python data are <0 or >1
    """
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
    """Tests that BBA values match between BioSNICAR data and the benchmark.

    Element-wise comparison between BBAs in equivalent positions in the Matlab benchmark
    dataset and the newly generated BioSNICAR dataset for the AD solver.

    Args:
        get_matlab_data: matlab-snicar-generated csv file of spectral and broadband albedo
        get_python_data: BioSNICAR generated csv file of spectral and broadband albedo
        set_tolerance: threshold error for BBAs to be considered equal

    Returns:
        None

    Raises:
        tests fail if the difference between any pair of BBA values exceeds set_tolerance

    """

    mat = get_matlab_data
    py = get_python_data
    tol = set_tolerance
    bb_py = py.loc[:, 481]
    bb_mat = mat.loc[:, 481]
    error = np.array(abs(bb_mat - bb_py))
    assert len(error[error > tol]) == 0


def test_compare_pyBBA_to_matBBA_clean(
    get_matlab_data_clean, get_python_data_clean, set_tolerance
):
    """Tests that BBA values match between BioSNICAR data and the benchmark for clean ice.

    Element-wise comparison between BBAs in equivalent positions in the Matlab benchmark
    dataset and the newly generated BioSNICAR dataset for the AD solver with no impurities.

    Args:
        get_matlab_data: matlab-snicar-generated csv file of spectral and broadband albedo
        get_python_data: BioSNICAR generated csv file of spectral and broadband albedo
        set_tolerance: threshold error for BBAs to be considered equal

    Returns:
        None

    Raises:
        tests fail if the difference between any pair of BBA values exceeds set_tolerance

    """

    mat = get_matlab_data_clean
    py = get_python_data_clean
    tol = set_tolerance
    bb_py = py.loc[:, 481]
    bb_mat = mat.loc[:, 481]
    error = np.array(abs(bb_mat - bb_py))
    assert len(error[error > tol]) == 0


def match_matlab_config(ice, illumination, rt_config, model_config, input_file):
    """Ensures model config is equal to the Matlab version used to generate benchmark data.

    This function resets values in instances of Ice, Illumination and ModelConfig to ensure
    equivalence between BioSNICAR and the Matlab code used to generate the benchmark data.
    Also ensures all vars have correct length, and re-executes the class functions in Ice and
    Illumination that update refractive indices and at-surface irradiance.

    Args:
        ice: instance of Ice class
        illumination: instance of Illumination class
        rt_config: instance of RTConfig class
        model_config: instance of ModelConfig class

    Returns:
        ice: updated instance of Ice class
        illumination: updated instance of Illumination class
        impurities: array of instances of Impurity class
        rt_config: updated instance of RTConfig class
        model_config: updated instance of ModelConfig class


    """

    nbr_lyr = 5
    # make sure ice config matches matlab benchmark
    ice.ri = 2
    ice.shp = [0] * nbr_lyr
    ice.shp_fctr = [0] * nbr_lyr
    ice.grain_ar = [0] * nbr_lyr
    ice.sfc = np.array([0.25] * model_config.nbr_wvl)
    ice.cdom = [0] * nbr_lyr
    ice.water = [0] * nbr_lyr
    ice.nbr_lyr = nbr_lyr
    ice.layer_type = [0] * nbr_lyr
    ice.rds = [ice.rds[0]] * nbr_lyr
    ice.rho = [ice.rho[0]] * nbr_lyr
    ice.lwc = [0] * nbr_lyr
    ice.lwc_pct_bbl = [0] * nbr_lyr
    ice.dz = [0.1] * nbr_lyr

    illumination.incoming = 4
    illumination.direct = 1

    # recalculate fluxes
    ice.calculate_refractive_index(input_file)
    illumination.calculate_irradiance()

    # make sure smoothing function is toggled off
    model_config.smooth = False

    # make sure impurities[0] is bc
    # (same bc used by matlab model)
    impurities = []

    conc = [0] * nbr_lyr
    impurity0 = Impurity("bc_ChCB_rn40_dns1270.npz", False, 0, "black_carbon", conc)
    impurities.append(impurity0)

    assert (impurities[0].name == "black_carbon") and (
        impurities[0].file == "bc_ChCB_rn40_dns1270.npz"
    )

    return ice, illumination, impurities, rt_config, model_config


def test_compare_pyspec_to_matspec_ad(get_matlab_data, get_python_data, set_tolerance):
    """Tests that spectral albedo values match between BioSNICAR data and the AD benchmark.

    Element-wise comparison between spectral albedo in equivalent positions in the Matlab benchmark
    dataset and the newly generated BioSNICAR dataset for the AD solver. Albedo is compared
    wavelength by wavelength for each column in the datasets.

    Args:
        get_matlab_data: matlab-snicar-generated csv file of spectral and broadband albedo
        get_python_data: BioSNICAR generated csv file of spectral and broadband albedo
        set_tolerance: threshold error for BBAs to be considered equal

    Returns:
        None

    Raises:
        tests fail if the difference between any pair of albedo values exceeds set_tolerance

    """

    mat = get_matlab_data
    py = get_python_data
    tol = set_tolerance
    py_spec = py.loc[:, :480]
    mat_spec = mat.loc[:, :480]
    error = np.array(abs(py_spec - mat_spec))
    assert len(error[error > tol]) == 0


def test_compare_pyspec_to_matspec_clean(
    get_matlab_data_clean, get_python_data_clean, set_tolerance
):
    """Tests that spectral albedo values match between BioSNICAR data and the Toon benchmark.

    Element-wise comparison between spectral albedo in equivalent positions in the Matlab benchmark
    dataset and the newly generated BioSNICAR dataset for the AD solver with no impurities.
    Albedo is compared wavelength by wavelength for each column in the datasets.

    Args:
        get_matlab_data_clean: matlab-snicar-generated csv file of spectral and broadband albedo
        get_python_data_clean: BioSNICAR generated csv file of spectral and broadband albedo
        set_tolerance: threshold error for BBAs to be considered equal

    Returns:
        None

    Raises:
        tests fail if the difference between any pair of albedo values exceeds set_tolerance

    """

    mat = get_matlab_data_clean
    py = get_python_data_clean
    tol = set_tolerance
    py_spec = py.loc[:, :480]
    mat_spec = mat.loc[:, :480]
    error = np.array(abs(py_spec - mat_spec))
    assert len(error[error > tol]) == 0


def test_plot_random_spectra_pairs(get_matlab_data, get_python_data, get_n_spectra):
    """Plots random selection of N spectra pairs and saves to /tests/test_data.

    Args:
        get_matlab_data: matlab-generated csv file of spectral and broadband albedo
        get_python_data: python-generated csv file of spectral and broadband albedo
        get_n_spectra: number of pairs of spectral albedo to plot

    Returns:
        None but saves py_mat_comparison.png to /tests/test_data/
    """

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


@pytest.mark.parametrize("dir", [0, 1])
@pytest.mark.parametrize("aprx", [1, 2, 3])
@pytest.mark.parametrize("inc", [0, 1, 2, 3, 4, 5, 6])
@pytest.mark.parametrize("ref", [0, 1, 2])
def test_config_fuzzer(dir, aprx, inc, ref, fuzz, input_file):
    """Checks model runs correctly with range of input value combinations.

    Fuzzer checks that model functions correctly across range of configurations.
    This fuzzer specifically checks rtm config and illumination parameters. The
    range of values is set in the parameterize decorators on this function and can
    be adjusted to test for specific failures or to increase test coverage. The
    defaults are designed to balance coverage with execution time. This func can be
    toggled off by setting fuzz to fase in conftest.py

    Args:
        dir: Boolean toggling between direct and diffuse irradiance
        aprx: choice of two-stream approximation
        inc: choice of spectral distribution of incoming irradiance
        ref: chocie of refractive indices (Warren 1984, Warren 2008, Picard 2016)
        fuzz: boolean toggling this fuxxing func on/off

    Returns:
        None

    Raises:
        tests fail if snicar throws an exception with a particular configuration
    """

    if fuzz:
        (
            ice,
            illumination,
            rt_config,
            model_config,
            plot_config,
            impurities,
        ) = setup_snicar("default")
        
        # make sure the BH83 data is used as per Matlab's version
        model_config.sphere_ice_path = "data/OP_data/480band/ice_spherical_grains_BH83/"
        model_config.bubbly_ice_path = "data/OP_data/480band/bubbly_ice_files_BH83/"

        ice, illumination, impurities, rt_config, model_config = match_matlab_config(
            ice, illumination, rt_config, model_config, input_file
        )

        rt_config.aprx_typ = aprx
        illumination.direct = dir
        ice.rf = ref
        ice.calculate_refractive_index(input_file)
        illumination.incoming = inc
        illumination.calculate_irradiance()

        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )
        outputs_toon = toon_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

        outputs_ad = adding_doubling_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config
        )

    else:
        pass

    return


@pytest.mark.parametrize("rds", [1000, 5000, 10000])
@pytest.mark.parametrize("rho", [400, 600, 800])
@pytest.mark.parametrize("zen", [50, 60, 70])
@pytest.mark.parametrize("dust", [0, 50000])
@pytest.mark.parametrize("algae", [0, 50000])
def test_var_fuzzer(rds, rho, zen, dust, algae, fuzz, input_file):
    """Checks model runs correctly with range of input value combinations.

    Fuzzer checks that model functions correctly across range of configurations.
    This fuzzer spoecifically checks input parameters including ice physical properties.
    The range of values is set in the parameterize decorators on this function and can
    be adjusted to test for specific failures or to increase test coverage. The
    defaults are designed to balance coverage with execution time. This func can be
    toggled off by setting fuzz to fase in conftest.py

    Args:
        rds: effective grain radius (um) of ice grains (lyr_typ==0) or air bubbles (lyr_typ==1)
        rho: density of ice layer (kg/m3)
        zen: zenith angle of direct beam (degrees from vertical)
        dust: concentration of mineral dust in each layer of the model (ppb)
        algae: concentration of glacir algae in each layer of the model (ppb)

    Returns:
        None

    Raises:
        tests fail if snicar throws an exception with a particular configuration
    """

    if fuzz:
        (
            ice,
            illumination,
            rt_config,
            model_config,
            plot_config,
            impurities,
        ) = setup_snicar("default")
        
        # make sure the BH83 data is used as per Matlab's version
        model_config.sphere_ice_path = "data/OP_data/480band/ice_spherical_grains_BH83/"
        model_config.bubbly_ice_path = "data/OP_data/480band/bubbly_ice_files_BH83/"

        ice, illumination, impurities, rt_config, model_config = match_matlab_config(
            ice, illumination, rt_config, model_config, input_file
        )

        impurities = []

        conc1 = [0] * len(ice.dz)
        conc1[0] = algae
        impurity0 = Impurity(
            "mie_sot_ChC90_dns_1317.npz",
            False,
            0,
            "black_carbon",
            conc1,
        )
        impurities.append(impurity0)

        conc2 = [0] * len(ice.dz)
        conc2[0] = dust
        impurity1 = Impurity(
            "dust_balkanski_central_size1.npz",
            False,
            0,
            "mineral_dust",
            conc2,
        )
        impurities.append(impurity1)

        ice.rds = [rds] * len(ice.dz)
        ice.rho = [rho] * len(ice.dz)
        illumination.solzen = zen
        illumination.calculate_irradiance()

        ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
        tau, ssa, g, L_snw = mix_in_impurities(
            ssa_snw, g_snw, mac_snw, ice, impurities, model_config
        )
        outputs_toon = toon_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
        )

        outputs_ad = adding_doubling_solver(
            tau, ssa, g, L_snw, ice, illumination, model_config
        )

    else:
        pass

    return


def test_tir_smoothing_high_sza():
    """Regression for #111: SG filter must not distort TIR albedo at SZA > ~55°.

    At oblique solar angles the direct beam can exceed the critical angle at
    ice anomalous-dispersion bands (~2.93–3.09 µm), producing physically exact
    albedo = 1.0 (total internal reflection, TIR).  The SG smoothing filter
    previously undershooting TIR bands (0.67–0.95 instead of 1.0) and falsely
    elevated adjacent guard bands (~0.33 instead of ~0.0004).

    Checks:
    - TIR bands remain exactly 1.0 after smoothing.
    - Guard bands (within window // 2 of TIR) deviate < 0.01 from raw values.
    - All albedo values lie in [0, 1].
    """
    from scipy.ndimage import binary_dilation

    input_file = "biosnicar/inputs.yaml"
    half_w = 7 // 2  # matches default window_size=7

    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar(
        input_file
    )
    # Solid ice with Fresnel surface — the only configuration that produces TIR
    ice.layer_type = [1]
    ice.dz = [0.02]
    ice.rds = [100]
    ice.rho = [300]
    ice.nbr_lyr = 1
    ice.lwc = [0]
    ice.lwc_pct_bbl = [0]
    ice.calculate_refractive_index(input_file)
    # Genuine TIR (n_re < 1 gate) fires where n_re < sin(SZA).  Min n_re for
    # Picard 2016 ice is ~0.954, so SZA must exceed arcsin(0.954) ≈ 72.5°.
    # Use SZA=80° to exercise the confirmed genuine-TIR path.
    illumination.solzen = 80
    illumination.calculate_irradiance()
    for imp in impurities:
        imp.conc = [0] * ice.nbr_lyr

    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(ssa_snw, g_snw, mac_snw, ice, impurities, model_config)

    model_config.smooth = False
    raw = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config).albedo.copy()

    model_config.smooth = True
    smoothed = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config).albedo.copy()

    tir = np.isclose(raw, 1.0)
    assert tir.sum() > 0, "Expected TIR bands at SZA=80 for solid ice (layer_type=1)"

    assert np.all(smoothed[tir] == 1.0), (
        f"TIR bands distorted by smoothing: min={smoothed[tir].min():.4f} "
        f"(expected 1.0; was as low as 0.67 before fix)"
    )

    guard = binary_dilation(tir, iterations=half_w) & ~tir
    if guard.any():
        max_dev = np.max(np.abs(smoothed[guard] - raw[guard]))
        assert max_dev < 0.01, (
            f"Guard bands falsely elevated: max deviation={max_dev:.4f} "
            f"(was ~0.33 before fix)"
        )

    assert np.all((smoothed >= 0.0) & (smoothed <= 1.0))


def test_tir_smoothing_no_regression_low_sza():
    """Regression for #111: at SZA=30 (no TIR), smoothing must be unchanged.

    The TIR guard logic must be a strict no-op when no TIR bands exist, so
    that standard snow configurations are unaffected by the fix.
    """
    from scipy.signal import savgol_filter

    input_file = "biosnicar/inputs.yaml"

    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar(
        input_file
    )
    illumination.solzen = 30
    illumination.calculate_irradiance()

    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(ssa_snw, g_snw, mac_snw, ice, impurities, model_config)

    model_config.smooth = False
    raw = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config).albedo.copy()

    model_config.smooth = True
    smoothed = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config).albedo.copy()

    assert not np.any(np.isclose(raw, 1.0)), "Unexpected TIR bands at SZA=30"

    # With no TIR bands the guard mask is empty, so output must equal plain SG + clip
    expected = np.clip(
        savgol_filter(raw, model_config.window_size, model_config.poly_order), 0.0, 1.0
    )
    np.testing.assert_array_equal(smoothed, expected)


def test_no_false_tir_sza89_4um():
    """Regression for false-TIR bug: no albedo=1 block at SZA=89, 4-5 µm.

    Before the fix, 70 bands at 4.145-4.835 µm (n_re≈1.34-1.35, n_im≈0.016-0.031)
    were overridden to Rf=1.0, creating a spurious second TIR block entirely
    outside the Reststrahlen region.  The physical Fresnel reflectance there is
    ~0.90 (grazing-angle high reflection), not 1.0.

    After the fix, those bands are computed correctly via the Fresnel formula and
    albedo stays near 0.90, not 1.0.
    """
    input_file = "biosnicar/inputs.yaml"
    ice, illumination, rt_config, model_config, plot_config, impurities = setup_snicar(
        input_file
    )
    ice.layer_type = [1]
    ice.dz = [0.02]
    ice.rds = [100]
    ice.rho = [300]
    ice.nbr_lyr = 1
    ice.lwc = [0]
    ice.lwc_pct_bbl = [0]
    ice.calculate_refractive_index(input_file)
    illumination.solzen = 89
    illumination.calculate_irradiance()
    for imp in impurities:
        imp.conc = [0] * ice.nbr_lyr

    ssa_snw, g_snw, mac_snw = get_layer_OPs(ice, model_config)
    tau, ssa, g, L_snw = mix_in_impurities(ssa_snw, g_snw, mac_snw, ice, impurities, model_config)

    model_config.smooth = False
    raw = adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config).albedo.copy()

    wvl = np.arange(0.205, 4.999, 0.01)
    false_tir_window = (wvl >= 4.1) & (wvl <= 4.9)

    # After fix: no band in 4.1-4.9 µm should be at 1.0 (n_re≈1.35 there)
    n_false = np.sum(np.isclose(raw[false_tir_window], 1.0))
    assert n_false == 0, (
        f"False TIR at SZA=89 in 4.1-4.9 µm: {n_false} bands incorrectly at 1.0 "
        "(n_re≈1.35 — outside Reststrahlen; see docs/TIR_CRITERION_BUG.md)"
    )
    # Bands should be near 0.85-0.95 (physically correct grazing-angle Fresnel)
    assert raw[false_tir_window].max() < 0.98, (
        f"Albedo unexpectedly high in 4.1-4.9 µm at SZA=89: {raw[false_tir_window].max():.4f}"
    )
    assert raw[false_tir_window].min() > 0.75, (
        f"Albedo unexpectedly low in 4.1-4.9 µm at SZA=89: {raw[false_tir_window].min():.4f}"
    )


if __name__ == "__main__":
    pass
