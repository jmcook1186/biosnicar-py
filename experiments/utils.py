"""
Joseph Cook, Aarhus University, Jan 2022

This script contains functions that a) generate 
SNICAR-predicted spectral albedo that approximate
field-measured spectra for a variety of weathering 
crust configurations; b) quantify the albedo change
resulting from a range of WC development scenarios 


includes:

find_best_params()
    the forward modelling script that retrieves the snicar params that generate the
    best-matching curve to a given field spectrum

call_snicar()
    the function used to do multiple snicar runs with params provided as a named tuple

match_field_spectra()
    the function for plotting simulated and measured field spectra in a multipanel fig

isolate_biological_effect()
    the function for calculating the albedo reduction due to bio vs phys processes by spectral differencing

build_LUT()
    function for constructing lookup table to be used in the inverse model

inverse_model()
    function finds the best matching entry in the LUTs for each field spectrum and returns the params
    used to generate it


"""
import collections as c
import sys

import dask
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr
from call_snicar import call_snicar

sys.path.append("./src")


def match_field_spectra(
    FIELD_DATA_FNAME,
    fnames,
    rho,
    rds,
    dz,
    alg,
    measured_cells,
    CIsites,
    LAsites,
    HAsites,
    APPLY_ARF,
    PLOT_ARF,
    ARF_CI,
    ARF_HA,
    SAVEPATH,
):

    """
    plot field against SNICAR spectra
    requires parameters to be known in advance and hard coded inside this function
    the relevant params can be generated using the find_best_params() func

    params:
    FIELD_DATA_FNAME = filename for spectral database

    The following params are parallel arrays - the order matters and
    must match the filenames! Pairs of values represent values for
    upper and lower layers in model. These are the values used to
    generate snicar spectrum to match field spectrum with name = fname[i]

    fnames = sample IDs for spectra to match model runs
    rho = pairs of density values [upper, lower]
    rds = pairs of r_eff values [upper, lower]
    dz = layer thickness [upper, lower] NB. upper is always 0.001
    alg = mass concentration of algae [upper, lower] NB. lower is always 0
    measured_cells = array of the actual measured cell concentration for each spectrum

    e.g.
    fnames= ['2016_WI_8','14_7_SB6','14_7_SB9','14_7_SB1','21_7_SB2','14_7_SB2', '22_7_SB3', 'RAIN']
    rho = [[550,550],[650,650],[800,800],[850,850],[750,750],[800,800],[800,800],[900,900]]
    rds = [[550,550],[650,650],[850,850],[850,850],[800,800],[800,800],[750,750],[900,900]]
    dz = [[0.001,0.3],[0.001,0.09],[0.001,0.03],[0.001,0.03],[0.001,0.02],[0.001,0.06],[0.001,0.05],[0.001,0.03]]
    alg = [[0,0],[0,0],[20000,0],[30000,0],[45000,0],[3000,0],[8000,0],[0,0]]

    returns:
    None, but saves figure to SAVEPATH

    """

    spectra = pd.read_csv(FIELD_DATA_FNAME)

    # reformat feld spectra to match snicar resolution
    spectra = spectra[::10]

    # gather spectra for each surface type
    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)]
    RAINspec = spectra["RAIN2"]

    if PLOT_ARF:
        plt.plot(
            spectra.Wavelength[0:100], ARF_CI[0:100], marker="x", label="clean ice ARF"
        ),
        plt.plot(
            spectra.Wavelength[0:100],
            ARF_HA[0:100],
            marker="o",
            linestyle="dashed",
            label="algal ice ARF",
        )
        plt.ylabel("Anitostropic Reflectance Factor"), plt.xlabel("Wavelength (nm)")
        plt.legend(loc="best")
        plt.savefig(str(SAVEPATH + "ARF.jpg"), dpi=300)

    # define local function for calling snicar
    def simulate_albedo(rds, rho, dz, alg):

        params = c.namedtuple(
            "params",
            "rho_layers, grain_rds, layer_type, dz,\
                 mss_cnc_glacier_algae, c_factor_GA, solzen",
        )
        params.grain_rds = rds
        params.rho_layers = rho
        params.layer_type = [1, 1]
        params.dz = dz
        params.mss_cnc_glacier_algae = alg
        params.c_factor_GA = 20
        params.solzen = 40
        albedo, BBA = call_snicar(params)

        return albedo, BBA

    # set up output array
    # and call snicar with each set of params
    OutArray = np.zeros(shape=(len(fnames), 480))

    for i in np.arange(0, len(fnames), 1):

        albedo, BBA = simulate_albedo(rds[i], rho[i], dz[i], alg[i])
        if APPLY_ARF:
            if alg[i][0] > 5000:
                albedo[15:230] = albedo[15:230] * ARF_HA
            else:
                albedo[15:230] = albedo[15:230] * ARF_CI

        OutArray[i, :] = albedo

    # calculate mean absolute error for model vs measured spectrum
    error = []
    for i in np.arange(0, len(fnames), 1):
        error.append(
            abs(np.mean(spectra[fnames[i]].iloc[0:130] - (OutArray[i, 15:145])))
        )

    # plot figure
    fig, axes = plt.subplots(4, 2, figsize=(10, 10))

    if APPLY_ARF:
        ylabel = "Reflectance"
        PlotName = "FieldvsMeasuredReflectance.jpg"
    else:
        ylabel = "Albedo"
        PlotName = "FieldvsMeasuredAlbedo.jpg"

    axes[0, 0].plot(spectra.Wavelength, spectra["23_7_SB1"], label="field")
    axes[0, 0].plot(
        spectra.Wavelength, OutArray[0, 15:230], label="model", linestyle="--"
    )
    axes[0, 0].set_ylim(0, 1), axes[0, 0].set_xlim(350, 1800)
    axes[0, 0].text(
        1400,
        0.2,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[0][0], rho[0][0], dz[0][1], error[0]
        ),
    )
    axes[0, 0].text(400, 0.1, "{} \n{} cells/mL ".format(fnames[0], measured_cells[0]))
    axes[0, 0].set_ylabel(ylabel), axes[0, 0].set_xlabel("Wavelength (nm)")
    axes[0, 0].legend(loc="best")

    axes[0, 1].plot(spectra.Wavelength, spectra["14_7_SB6"])
    axes[0, 1].plot(spectra.Wavelength, OutArray[1, 15:230], linestyle="--")
    axes[0, 1].set_ylim(0, 1), axes[0, 1].set_xlim(350, 1800)
    axes[0, 1].text(
        1450,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[1][0], rho[1][0], dz[1][1], error[1]
        ),
    )
    axes[0, 1].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[1], measured_cells[1]))
    axes[0, 1].set_ylabel(ylabel), axes[0, 1].set_xlabel("Wavelength (nm)")

    axes[1, 0].plot(spectra.Wavelength, spectra["14_7_SB9"])
    axes[1, 0].plot(spectra.Wavelength, OutArray[2, 15:230], linestyle="--")
    axes[1, 0].set_ylim(0, 1), axes[1, 0].set_xlim(350, 1800)
    axes[1, 0].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[2][0], rho[2][0], dz[2][1], error[2]
        ),
    )
    axes[1, 0].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[2], measured_cells[2]))
    axes[1, 0].set_ylabel(ylabel), axes[1, 0].set_xlabel("Wavelength (nm)")

    axes[1, 1].plot(spectra.Wavelength, spectra["14_7_SB1"])
    axes[1, 1].plot(spectra.Wavelength, OutArray[3, 15:230], linestyle="--")
    axes[1, 1].set_ylim(0, 1), axes[1, 1].set_xlim(350, 1800)
    axes[1, 1].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[3][0], rho[3][0], dz[3][1], error[3]
        ),
    )
    axes[1, 1].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[3], measured_cells[3]))
    axes[1, 1].set_ylabel(ylabel), axes[1, 1].set_xlabel("Wavelength (nm)")

    axes[2, 0].plot(spectra.Wavelength, spectra["22_7_SB5"])
    axes[2, 0].plot(spectra.Wavelength, OutArray[4, 15:230], linestyle="--")
    axes[2, 0].set_ylim(0, 1), axes[2, 0].set_xlim(350, 1800)
    axes[2, 0].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[4][0], rho[4][0], dz[4][1], error[4]
        ),
    )
    axes[2, 0].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[4], measured_cells[4]))
    axes[2, 0].set_ylabel(ylabel), axes[2, 0].set_xlabel("Wavelength (nm)")

    axes[2, 1].plot(spectra.Wavelength, spectra["14_7_SB2"])
    axes[2, 1].plot(spectra.Wavelength, OutArray[5, 15:230], linestyle="--")
    axes[2, 1].set_ylim(0, 1), axes[2, 1].set_xlim(350, 1800)
    axes[2, 1].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[5][0], rho[5][0], dz[5][1], error[5]
        ),
    )
    axes[2, 1].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[5], measured_cells[5]))
    axes[2, 1].set_ylabel(ylabel), axes[2, 1].set_xlabel("Wavelength (nm)")

    axes[3, 0].plot(spectra.Wavelength, spectra["22_7_SB3"])
    axes[3, 0].plot(spectra.Wavelength, OutArray[6, 15:230], linestyle="--")
    axes[3, 0].set_ylim(0, 1), axes[3, 0].set_xlim(350, 1800)
    axes[3, 0].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[6][0], rho[6][0], dz[6][1], error[6]
        ),
    )
    axes[3, 0].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[6], measured_cells[6]))
    axes[3, 0].set_ylabel(ylabel), axes[3, 0].set_xlabel("Wavelength (nm)")

    axes[3, 1].plot(spectra.Wavelength, RAINspec)
    axes[3, 1].plot(spectra.Wavelength, OutArray[7, 15:230], linestyle="--")
    axes[3, 1].set_ylim(0, 1), axes[3, 1].set_xlim(350, 1800)
    axes[3, 1].text(
        1400,
        0.3,
        "r_eff: {}\nrho: {}\ndz: {}\nC-factor: 30\nerror: {:.3f}".format(
            rds[7][0], rho[7][0], dz[7][1], error[7]
        ),
    )
    axes[3, 1].text(400, 0.8, "{}: \n{} cells/mL".format(fnames[1], measured_cells[1]))
    axes[3, 1].set_ylabel(ylabel), axes[3, 1].set_xlabel("Wavelength (nm)")

    fig.tight_layout()
    plt.savefig(str(SAVEPATH + PlotName), dpi=300)

    return True


def build_LUT(
    solzen,
    dz,
    densities,
    radii,
    algae,
    wavelengths,
    save_LUT,
    APPLY_ARF,
    ARF_CI,
    ARF_HA,
    SAVEPATH,
):

    """
    generates LUTs used to invert BioSNICAR in RISA project

    params:
    ice_rds: fixed effective bubble radius for solid ice layers (default = 525)
    ice dens: fixed density for solid ice layers (default = 894)
    zeniths: range of solar zenith angles to loop over
    dz: thickness of each vertical layer
    densities: densities for top layer. Lower layers predicted by exponential model
    algae: mass mixing ratio of algae in top layer
    wavelengths: wavelength range, default is np.arange(0.2, 5, 0.01)
    save_LUT: Boolean to toggle saving to npy file
    SAVEPATH: directory to save LUT

    returns:
    WCthickLUT: for each index position in the spectraLUT, this holds the WC
                thickness in the corresponding index position
    SpectraLUT: ND array containing 480element spectrum for each
                dens/alg/zen combination

    return spectraLUT


    """
    LUT = []

    @dask.delayed
    def run_sims(dens, rad, dz, alg, zen):

        params = c.namedtuple(
            "params",
            "rho_layers, grain_rds, layer_type, c_factor_GA, dz, mss_cnc_glacier_algae, solzen",
        )
        params.rho_layers = [dens, dens]
        params.grain_rds = [rad, rad]  # set equal to density
        params.layer_type = [1, 1]
        params.dz = [0.001, dz]
        params.mss_cnc_glacier_algae = [alg, 0]
        params.solzen = zen
        params.c_factor_GA = 20

        albedo, BBA = call_snicar(params)

        return albedo

    for z in np.arange(0, len(solzen), 1):
        for i in np.arange(0, len(densities), 1):
            for j in np.arange(0, len(radii), 1):
                for p in np.arange(0, len(dz), 1):
                    for q in np.arange(0, len(algae), 1):

                        albedo = run_sims(
                            densities[i], radii[j], dz[p], algae[q], solzen[z]
                        )

                        LUT.append(albedo)

    LUT = dask.compute(*LUT, num_workers=12)
    LUT = np.array(LUT).reshape(
        len(solzen), len(densities), len(radii), len(dz), len(algae), len(wavelengths)
    )

    # move the ARF application to new loop because dask compute objets are immutable
    # i.e. modifications to albedo must be done post-compute
    if APPLY_ARF:
        for z in np.arange(0, len(solzen), 1):
            for i in np.arange(0, len(densities), 1):
                for j in np.arange(0, len(radii), 1):
                    for p in np.arange(0, len(dz), 1):
                        for q in np.arange(0, len(algae), 1):

                            if algae[q] > 5000:

                                LUT[z, i, j, p, q, 0:215] = (
                                    LUT[z, i, j, p, q, 0:215] * ARF_HA
                                )

                            else:
                                LUT[z, i, j, p, q, 0:215] = (
                                    LUT[z, i, j, p, q, 0:215] * ARF_CI
                                )

    if save_LUT:
        np.save(str(SAVEPATH + "LUT.npy"), LUT)

    return LUT


def inverse_model(
    FIELD_DATA_FNAME, LUT_PATH, DZ, DENSITIES, RADII, ZENS, ALGAE, WAVELENGTHS
):

    # read in and reshape luts and field spectra
    field_data = pd.read_csv(FIELD_DATA_FNAME, index_col=None)
    field_spectra = field_data[::10]
    lut = np.load(str(LUT_PATH + "LUT.npy"))
    vis_start_idx = 15
    vis_end_idx = 55
    nir_start_idx = 55
    nir_end_idx = 230
    lut_nir = lut[:, :, :, :, :, nir_start_idx:nir_end_idx]
    lut_vis = lut[:, :, :, :, :, vis_start_idx:vis_end_idx]

    #   LUT = np.array(LUT).reshape(
    #         len(solzen),
    #   len(densities),
    #   len(radii),
    #   len(dz),
    #   len(algae),
    #   len(wavelengths)
    #     )
    # TODO: make this dynamic depending on var
    # values from config
    flat_nir_lut = lut_nir.reshape(
        len(ZENS) * len(DENSITIES) * len(RADII) * len(DZ) * len(ALGAE), 175
    )

    output = pd.DataFrame()

    retrieved_zen = []
    retrieved_density = []
    retrieved_radii = []
    retrieved_dz = []
    retrieved_algae = []
    names = []
    errors = []

    for (name, data) in field_spectra.iteritems():

        names.append(name)

        # step 1: find params that match best in NIR
        error_array = abs(flat_nir_lut - np.array(data[40:]))
        error_mean = np.mean(error_array, axis=1)
        index = np.argmin(error_mean)
        param_idx_phys = np.unravel_index(index, [len(ZENS) * len(DENSITIES) * len(RADII) * len(DZ) * len(ALGAE), 1])

        # step 2: fix physical params and minimise error in visd by varying algae
        lut2 = lut_vis[
            param_idx_phys[0],
            param_idx_phys[1],
            param_idx_phys[2],
            param_idx_phys[3],
            :,
            :,
        ]
        error_array = abs(lut2 - np.array(data[vis_start_idx:vis_end_idx]))
        error_mean = np.mean(error_array, axis=1)
        index = np.argmin(error_mean)
        param_idx_alg = np.unravel_index(index, [1, 1, 1, 1, len(ALGAE), 1])

        # setp 3:
        retrieved_zen.append(ZENS[param_idx_phys[0]])
        retrieved_density.append(DENSITIES[param_idx_phys[1]])
        retrieved_radii.append(RADII[param_idx_phys[2]])
        retrieved_dz.append(DZ[param_idx_phys[3]])
        retrieved_algae.append(ALGAE[param_idx_alg[4]])

        total_error = abs(
            lut[
                param_idx_phys[0],
                param_idx_phys[1],
                param_idx_phys[2],
                param_idx_phys[3],
                param_idx_alg[4],
                vis_start_idx:nir_end_idx,
            ]
            - data
        )
        errors.append(np.mean(total_error))

    output["fname"] = names
    output["solzen"] = retrieved_zen
    output["density"] = retrieved_density
    output["radii"] = retrieved_radii
    output["dz"] = retrieved_dz
    output["algae"] = retrieved_algae
    output["error"] = errors

    return output


def isolate_biological_effect(FIELD_DATA_FNAME, CIsites, LAsites, HAsites, SAVEPATH):

    """
    This function estimates the albedo reduction resulting from the ice physical changes
    versus the biological growth.

    Some nuance to the interpretation because the ce surface likely would not
    degrade to the same extent without the algal bloom.

    """

    # read in spectral database
    spectra = pd.read_csv(FIELD_DATA_FNAME)

    # reformat feld spectra to match snicar resolution
    spectra = spectra[::10]
    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)]
    meanCI = CIspec.mean(axis=1)
    meanHA = HAspec.mean(axis=1)
    meanLA = LAspec.mean(axis=1)

    # define local function for calling snicar
    def simulate_albedo(rds, rho, dz, alg):
        params = c.namedtuple(
            "params",
            "rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, solzen",
        )
        params.grain_rds = rds
        params.rho_layers = rho
        params.layer_type = [1, 1]
        params.dz = dz
        params.mss_cnc_glacier_algae = alg
        params.solzen = 53
        albedo, BBA = call_snicar(params)
        return albedo, BBA

    # call snicar to generat esimulated spectrum for LA and HA using predetermined params

    SNICARalbedoLA, BBA = simulate_albedo([800, 800], [800, 800], [0.001, 0.1], [0, 0])
    SNICARalbedoHA, BBA = simulate_albedo([900, 900], [800, 800], [0.001, 0.08], [0, 0])
    SNICARalbedoLA = SNICARalbedoLA[15:230]
    SNICARalbedoHA = SNICARalbedoHA[15:230]

    # plot figure
    x = np.arange(350, 2500, 10)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.plot(x, meanCI, linestyle="--", alpha=0.4, label="Clean ice (mean)")
    ax1.plot(x, meanLA, linestyle="-.", alpha=0.4, label="Algal ice (mean)")
    ax1.plot(
        x, SNICARalbedoLA, linestyle="dotted", alpha=0.4, label="Clean ice (model)"
    )
    ax1.fill_between(x, meanCI, SNICARalbedoLA, alpha=0.2)
    ax1.fill_between(x, SNICARalbedoLA, meanLA, color="k", alpha=0.2)
    ax1.set_xlim(350, 1500), ax1.legend(loc="best")
    ax1.set_ylabel("Albedo"), ax1.set_xlabel("Wavelength (nm)")

    ax2.plot(x, meanCI, linestyle="--", alpha=0.4, label="Clean ice (mean)")
    ax2.plot(x, meanHA, linestyle="-.", alpha=0.4, label="Algal ice (mean)")
    ax2.plot(
        x, SNICARalbedoHA, linestyle="dotted", alpha=0.4, label="Clean ice (model)"
    )
    ax2.fill_between(x, meanCI, SNICARalbedoHA, alpha=0.2)
    ax2.fill_between(x, SNICARalbedoHA, meanHA, color="k", alpha=0.2)
    ax2.set_xlim(350, 1500), ax2.legend(loc="best")
    ax2.set_ylabel("Albedo"), ax2.set_xlabel("Wavelength (nm)")
    fig.tight_layout()
    plt.savefig(str(SAVEPATH + "/BiovsPhysEffect.jpg"), dpi=300)

    # define incoming to calculate broadband albedo
    incoming = xr.open_dataset(
        "/home/joe/Code/BioSNICAR_GO_PY/Data/Mie_files/480band/fsds/swnb_480bnd_sas_clr_SZA60.nc"
    )
    incoming = incoming["flx_frc_sfc"].values
    incoming = incoming[15:230]

    # calculate broadband albedo of each case
    LA_BBA = np.sum(meanLA * incoming) / np.sum(incoming)
    HA_BBA = np.sum(meanHA * incoming) / np.sum(incoming)
    CI_BBA = np.sum(meanCI * incoming) / np.sum(incoming)
    CI2_BBA_LA = np.sum(SNICARalbedoLA * incoming) / np.sum(incoming)
    CI2_BBA_HA = np.sum(SNICARalbedoHA * incoming) / np.sum(incoming)

    # calculate change due to bio/phys as BBA difference
    delAbioLA = CI2_BBA_LA - LA_BBA
    delAphysLA = CI_BBA - CI2_BBA_LA
    delAbioHA = CI2_BBA_HA - HA_BBA
    delAphysHA = CI_BBA - CI2_BBA_HA

    return delAbioLA, delAphysLA, delAbioHA, delAphysHA


def run_best_params(
    SAVEPATH,
    ALL_FIELD_SAMPLES,
    FIELD_DATA_FNAME,
    CI_SITES,
    LA_SITES,
    HA_SITES,
    WEIGHT,
    CLEAN,
):
    """
    function calls out to find_best_params
    """
    ResultArray = np.zeros(shape=(len(ALL_FIELD_SAMPLES), 6))

    for i in np.arange(0, len(ALL_FIELD_SAMPLES), 1):

        fn = ALL_FIELD_SAMPLES[i]
        Results = find_best_params(
            FIELD_DATA_FNAME, fn, CI_SITES, LA_SITES, HA_SITES, WEIGHT, CLEAN
        )
        Results = np.array(Results)
        best_idx = Results[:, 1].argmin()
        best_params = Results[best_idx, 0]
        best_error = Results[best_idx, 1]
        best_dens, best_rds, best_dz, best_alg, best_zen = zip(best_params)
        ResultArray[i, :] = (
            np.array(best_dens),
            np.array(best_rds),
            np.array(best_dz),
            np.array(best_alg),
            np.array(best_zen),
            np.array(best_error),
        )

    Out = pd.DataFrame(
        data=ResultArray, columns=["dens", "rds", "dz", "alg", "zen", "spec_err"]
    )
    Out.index = ALL_FIELD_SAMPLES
    Out.to_csv(str(SAVEPATH + "retrieved_params.csv"))

    return True


def find_best_params(
    FIELD_DATA_FNAME, sampleID, CIsites, LAsites, HAsites, weight, clean=True
):

    """
    This function will return the SNICAR parameter set that provides the
    best approximation to a given field spectrum.

    """

    spectra = pd.read_csv(FIELD_DATA_FNAME)
    spectra = spectra[::10]

    CIspec = spectra[spectra.columns.intersection(CIsites)]
    HAspec = spectra[spectra.columns.intersection(HAsites)]
    LAspec = spectra[spectra.columns.intersection(LAsites)]

    if sampleID == "CImean":
        field_spectrum = CIspec.mean(axis=1)
    elif sampleID == "HAmean":
        field_spectrum = HAspec.mean(axis=1)
    elif sampleID == "LAmean":
        field_spectrum = LAspec.mean(axis=1)
    elif sampleID == "RAIN":
        field_spectrum = spectra["RAIN2"]
    else:
        field_spectrum = spectra[sampleID]

    # calculate 2BDA index of field spectrum
    BDA2idx = np.array(field_spectrum)[36] / np.array(field_spectrum)[31]

    dens = [550, 600, 650, 700, 750, 800, 850]
    dz = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2]
    rds = [600, 700, 800, 900, 1000]
    alg = [0, 2500, 5000, 7500, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000]
    c_factors = [10, 20, 30]
    solzen = [40, 45, 50]

    @dask.delayed
    def run_sims(i, j, k, p, q, z):

        params = c.namedtuple(
            "params",
            "rho_layers, grain_rds, layer_type, dz, mss_cnc_glacier_algae, c_factor_GA, solzen",
        )
        params.rho_layers = [i, i]
        params.grain_rds = [j, j]
        params.layer_type = [1, 1]
        params.dz = [0.001, k]
        params.mss_cnc_glacier_algae = [p, 0]
        params.c_factor_GA = c_factors[q]
        params.solzen = z

        albedo, BBA = call_snicar(params)

        error_vis = np.mean(abs(albedo[15:55] - field_spectrum[0:40]))
        error_nir = np.mean(abs(albedo[55:100] - field_spectrum[40:85]))
        error = ((error_vis * weight) + error_nir) / (1 + weight)

        params = (i, j, k, p, z)
        out = (params, error)

        return out

    Out = []
    # now use the reduced LUT to call snicar and obtain best matching spectrum
    for i in dens:
        for j in rds:
            for k in dz:
                for p in alg:
                    for z in solzen:
                        out = run_sims(i, j, k, p, z)
                        Out.append(out)

    Result = dask.compute(*Out, num_workers=12)

    return Result


def BDA2_of_field_samples():

    """
    2DBA index calculated from field samples after field spectra are averaged over S2 band 4 and 5 wavelengths
    weighted by the sensor spectral response function for each band. The index is then calculated as B5/B4
    and the cell concentration predicted using Wang et al's (2018) conversion equation.

    """

    spectra = pd.read_csv(
        "/home/joe/Code/Remote_Ice_Surface_Analyser/Training_Data/HCRF_master_16171819.csv"
    )

    # reformat LUT: flatten LUT from 3D to 2D array with one column per combination
    # of RT params, one row per wavelength

    responsefunc = pd.read_csv(
        "/home/joe/Code/Remote_Ice_Surface_Analyser/S2SpectralResponse.csv"
    )
    func04 = responsefunc["B4"].loc[
        (responsefunc["SR_WL"] > 650) & (responsefunc["SR_WL"] <= 680)
    ]
    func05 = responsefunc["B5"].loc[
        (responsefunc["SR_WL"] > 698) & (responsefunc["SR_WL"] <= 713)
    ]

    filenames = []
    Idx2DBAList = []
    prd2DBAList = []
    Idx2DBA_S2List = []
    prd2DBA_S2List = []
    Idx2DBA_Ideal_List = []
    prd2DBA_Ideal_List = []

    for i in np.arange(0, len(spectra.columns), 1):

        if i != "Wavelength":

            colname = spectra.columns[i]
            spectrum = np.array(spectra[colname])

            B04 = np.mean(spectrum[300:330] * func04)
            B05 = np.mean(spectrum[348:363] * func05)
            Idx2DBA = spectrum[355] / spectrum[315]
            prd2DBA = 10e-35 * Idx2DBA * np.exp(87.015 * Idx2DBA)
            Idx2DBA_S2 = B05 / B04
            prd2DBA_S2 = 10e-35 * Idx2DBA_S2 * np.exp(87.015 * Idx2DBA_S2)
            Idx2DBA_Ideal = spectrum[360] / spectrum[330]
            prd2DBA_Ideal = 10e-35 * Idx2DBA_Ideal * np.exp(87.015 * Idx2DBA_Ideal)

            filenames.append(colname)
            Idx2DBAList.append(Idx2DBA)
            prd2DBAList.append(prd2DBA)
            Idx2DBA_S2List.append(Idx2DBA_S2)
            prd2DBA_S2List.append(prd2DBA_S2)
            Idx2DBA_Ideal_List.append(Idx2DBA_Ideal)
            prd2DBA_Ideal_List.append(prd2DBA_Ideal)

    Out = pd.DataFrame()
    Out["filename"] = filenames
    Out["2DBAIdx"] = Idx2DBAList
    Out["2DBAPrediction"] = prd2DBAList
    Out["2DBA_S2Idx"] = Idx2DBA_S2List
    Out["2DBA_S2Prediction"] = prd2DBA_S2List
    Out["2DBAIdx_Ideal"] = Idx2DBA_Ideal_List
    Out["2DBAPrediction_Ideal"] = prd2DBA_Ideal_List

    return Out


def compare_predicted_and_measured(SAVEPATH, path_to_metadata):

    ## imports and data organisation
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    DF = pd.read_csv(path_to_metadata)

    measured_cells = DF["measured_cells"]
    modelled_cells = DF["algae_cells_inv_model_S2"]
    BDA2_cells = DF["cells_BDA2_centre_wang"]
    BDA2_idx = DF["BDA2_centre"]

    ## regression models

    # Ordinary least squares regression
    model1 = sm.OLS(modelled_cells, measured_cells).fit()
    summary1 = model1.summary()
    test_x = [
        0,
        1000,
        5000,
        7500,
        10000,
        12500,
        15000,
        17500,
        20000,
        25000,
        30000,
        35000,
        40000,
        50000,
    ]
    ypred1 = model1.predict(test_x)

    # regress measured cells against band index
    # use this to give predictive linear model

    BDA2_PredModel = sm.OLS(measured_cells, sm.add_constant(BDA2_idx)).fit()
    BDA2_PredModel_r2 = np.round(BDA2_PredModel.rsquared, 3)
    BDA2_PredModel_y = BDA2_PredModel.predict(sm.add_constant(BDA2_idx))

    # regress BDA2 predicted cells against measured cells
    model2 = sm.OLS(BDA2_PredModel_y, measured_cells).fit()
    summary2 = model2.summary()
    test_x = [
        0,
        1000,
        5000,
        7500,
        10000,
        12500,
        15000,
        17500,
        20000,
        25000,
        30000,
        35000,
        40000,
        50000,
    ]
    ypred2 = model2.predict(test_x)

    # multipanel figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    ax1.plot(
        measured_cells, color="k", marker="x", linestyle="None", label="field-measured"
    )
    ax1.plot(
        modelled_cells,
        color="b",
        marker="o",
        markerfacecolor="None",
        alpha=0.6,
        linestyle="None",
        label="RTM model prediction",
    )
    ax1.plot(
        BDA2_PredModel_y,
        color="r",
        marker="^",
        markerfacecolor="r",
        alpha=0.3,
        linestyle="None",
        label="new 2BDA model prediction",
    )
    ax1.set_ylabel("Algal concentration (cells/mL)")
    ax1.set_xticks(range(len(measured_cells)))
    ax1.set_xticklabels([])
    ax1.legend(loc="upper left")
    ax1.set_xlabel("Individual samples")
    ax1.set_ylim(0, 65000)

    ax2.scatter(
        measured_cells,
        modelled_cells,
        marker="o",
        facecolor="None",
        color="b",
        alpha=0.6,
        label="RTM\nr$^2$ = {}\np = {}".format(
            np.round(model1.rsquared, 3), np.round(model1.pvalues[0], 8)
        ),
    )
    ax2.plot(test_x, ypred1, linestyle="dotted", color="b", alpha=0.6)
    ax2.scatter(
        measured_cells,
        BDA2_PredModel_y,
        marker="^",
        facecolor="r",
        color="r",
        alpha=0.3,
        label="2BDA\nr$^2$ = {}\np = {}".format(
            np.round(model2.rsquared, 3), np.round(model2.pvalues[0], 10)
        ),
    )
    ax2.plot(test_x, ypred2, linestyle="dashed", color="r", alpha=0.3)
    ax2.set_ylabel("Algal concentration,\n cells/mL (field)")
    ax2.set_xlabel("Algal concentration,\n clls/mL (predicted by model)")
    ax2.set_xlim(0, 50000), ax2.set_ylim(0, 60000)

    ax2.legend(loc="upper left")

    fig.tight_layout()

    SAVEPATH = "/home/joe/Code/Remote_Ice_Surface_Analyser/Manuscript/Figures"
    fig.savefig(str(SAVEPATH + "/measured_modelled_algae.png"), dpi=300)

    return
