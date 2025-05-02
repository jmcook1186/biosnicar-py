#!/usr/bin/python
"""Radiative transfer solver using the adding-doubling method.

Here the ice/impurity optical properties and illumination conditions
are used to calculate energy fluxes between the ice, atmosphere and
underlying substrate.

Typically this function would be called from snicar_driver() because
it takes as inputs intermediates that are calculated elsewhere.
Specifically, the functions setup_snicar(), get_layer_OPs() and
mix_in_impurities() are called to generate tau, ssa, g and L_snw,
which are then passed as inputs to adding_doubling_solver().

The adding-doubling routine implemented here originates in Brieglib
and Light (2007) and was coded up in Matlab by Chloe Whicker and Mark
Flanner to be published in Whicker (2022: The Cryosphere). Their
scripts were the jumping off point for this script and their code is still
used to benchmark this script against.

The adding-doubling solver implemented here has been shown to do an excellent
job at simulating solid glacier ice. This solver can either treat ice as a
"granular" material with a bulk medium of air with discrete ice grains, or
as a bulk medium of ice with air inclusions. In the latter case, the upper
boundary is a Fresnel reflecting surface. Total internal reflection is accounted
for if the irradiance angles exceed the critical angle.

This is always the appropriate solver to use in any  model configuration where
solid ice layers and fresnel reflection are included.

"""

import numpy as np
from scipy.signal import savgol_filter

from biosnicar.classes.outputs import Outputs


def adding_doubling_solver(tau, ssa, g, L_snw, ice, illumination, model_config):
    """control function for the adding-doubling solver.

    Makes function calls in sequence to generate, then return, an instance of
    Outputs class.

    Args:
        tau: optical thickness of ice column in m/m
        ssa: single scattering albedo of ice column (dimensionless)
        g: asymmetry parameter for ice column
        L_snw: mass of ice in lg
        ice: instance of Ice class
        illumination: instance of Illumination class
        model_config: instance of ModelConfig class

    Returns:
        outputs: Instance of Outputs class

    Raises:
        ValueError if violation of conservation of energy detected

    """

    # DEFINE CONSTANTS AND ARRAYS
    (
        tau0,
        g0,
        ssa0,
        epsilon,
        exp_min,
        nr,
        mu0,
        mu0n,
        trnlay,
        rdif_a,
        rdif_b,
        tdif_a,
        tdif_b,
        rdir,
        tdir,
        lyrfrsnl,
        trnlay,
        rdif_a,
        rdif_b,
        tdif_a,
        tdif_b,
        rdir,
        tdir,
        rdndif,
        trntdr,
        trndir,
        trndif,
        fdirup,
        fdifup,
        fdirdn,
        fdifdn,
        dfdir,
        dfdif,
        F_up,
        F_dwn,
        F_abs,
        F_abs_vis,
        F_abs_nir,
        rupdif,
        rupdir,
    ) = define_constants_arrays(tau, g, ssa, illumination, ice, model_config)

    # proceed down one layer at a time: if the total transmission to
    # the interface just above a given layer is less than trmin, then no
    # Delta-Eddington computation for that layer is done.

    for lyr in np.arange(0, ice.nbr_lyr, 1):  # loop through layers
        # condition: if current layer is above fresnel layer or the
        # top layer is a Fresnel layer
        if lyr < lyrfrsnl:
            mu0n = mu0

        else:
            # within or below fl
            mu0n = mu0n

        rdir, tdir, ts, ws, gs, lm = calc_reflectivity_transmittivity(
            tau0,
            ssa0,
            g0,
            lyr,
            model_config,
            exp_min,
            rdif_a,
            tdif_a,
            trnlay,
            mu0n,
            epsilon,
            rdir,
            tdir,
        )
        # recalculate rdif,tdif using direct angular integration over rdir,tdir,
        # since Delta-Eddington rdif formula is not well-behaved (it is usually
        # biased low and can even be negative)  use ngmax angles and gaussian
        # integration for most accuracy:
        smt, smr, swt = apply_gaussian_integral(
            model_config, exp_min, ts, ws, gs, epsilon, lm, lyr, rdif_a, tdif_a
        )

        rdif_a, tdif_a, rdif_b, tdif_b = update_transmittivity_reflectivity(
            swt, smr, smt, lyr, rdif_a, tdif_a, rdif_b, tdif_b
        )

        if lyr == lyrfrsnl:
            (
                rdif_a,
                rdif_b,
                tdif_a,
                tdif_b,
                trnlay,
                rdir,
                tdir,
            ) = calc_correction_fresnel_layer(
                model_config,
                ice,
                illumination,
                mu0n,
                mu0,
                nr,
                rdif_a,
                rdif_b,
                tdif_a,
                tdif_b,
                trnlay,
                lyr,
                rdir,
                tdir,
            )

        trndir, trntdr, rdndif, trndif = calc_reflection_transmission_from_top(
            lyr,
            trnlay,
            rdif_a,
            rdir,
            tdif_a,
            rdif_b,
            tdir,
            tdif_b,
            model_config,
            ice,
            trndir,
            rdndif,
            trntdr,
            trndif,
        )

    rupdif, rupdir = calc_reflection_below(
        ice,
        model_config,
        rdif_a,
        rdif_b,
        tdif_a,
        tdif_b,
        trnlay,
        rdir,
        tdir,
        rupdif,
        rupdir,
    )

    fdirup, fdifup, fdirdn, fdifdn = trans_refl_at_interfaces(
        model_config,
        ice,
        rupdif,
        rupdir,
        rdndif,
        trndir,
        trndif,
        trntdr,
        fdirup,
        fdirdn,
        fdifup,
        fdifdn,
        dfdir,
        dfdif,
    )

    albedo, F_abs, F_btm_net, F_top_pls = calculate_fluxes(
        model_config,
        ice,
        illumination,
        fdirup,
        fdifup,
        fdirdn,
        fdifdn,
        F_up,
        F_dwn,
        F_abs,
        F_abs_vis,
        F_abs_nir,
    )

    conservation_of_energy_check(illumination, F_abs, F_btm_net, F_top_pls)

    outputs = get_outputs(illumination, albedo, model_config, L_snw, F_abs, F_btm_net)

    if model_config.smooth:
        outputs.albedo = apply_smoothing_function(outputs.albedo, model_config)

    return outputs


def calc_reflectivity_transmittivity(
    tau0,
    ssa0,
    g0,
    lyr,
    model_config,
    exp_min,
    rdif_a,
    tdif_a,
    trnlay,
    mu0n,
    epsilon,
    rdir,
    tdir,
):
    """Calculates reflectivity and transmissivity.

    Sets up new variables, applies delta transformation and makes
    initial calculations of reflectivity and transmissivity in each
    layer.

    Args:
        tau0: initial optical thickness
        ssa0: initial single scattering albedo
        g0: initial asymmetry parameter
        lyr: index of current layer
        model_config: instance of ModelConfig class
        exp_min: small number to avoid /zero error
        rdif_a: reflectivity to diffuse irradiance at polarization angle == perpendicular
        tdif_a: transmissivity to diffuse irradiance at polarization angle == perpendicular
        trnlay: transmission through layer
        mu0n: incident beam angle adjusted for refraction
        epsilon: small number to avoid singularity
        rdir: reflectivity to direct beam
        tdir: transmissivity to direct beam

    Returns:
        rdir:
        tdir:
        ts:
        ws:
        gs:
        lm:

    """
    # calculation over layers with penetrating radiation
    # includes optical thickness, single scattering albedo,
    # asymmetry parameter and total flux
    tautot = tau0[:, lyr]
    wtot = ssa0[:, lyr]
    gtot = g0[:, lyr]
    ftot = g0[:, lyr] * g0[:, lyr]

    # coefficient for delta eddington solution for all layers
    # Eq. 50: Briegleb and Light 2007
    # layer delta-scaled extinction optical depth
    ts = (1 - (wtot * ftot)) * tautot
    ws = ((1 - ftot) * wtot) / (
        1 - (wtot * ftot)
    )  # layer delta-scaled single scattering albedo
    gs = gtot / (1 + gtot)  # layer delta-scaled asymmetry parameter
    lm = np.sqrt(3 * (1 - ws) * (1 - ws * gs))  # lambda
    ue = (
        1.5 * (1 - ws * gs) / lm
    )  # u equation, term in diffuse reflectivity and transmissivity
    extins = np.maximum(
        np.full((model_config.nbr_wvl,), exp_min), np.exp(-lm * ts)
    )  # extinction, MAX function lyr keeps from getting an error
    # if the exp(-lm*ts) is < 1e-5
    ne = (ue + 1) ** 2 / extins - (
        ue - 1
    ) ** 2 * extins  # N equation, term in diffuse reflectivity and transmissivity

    # ! first calculation of rdif, tdif using Delta-Eddington formulas
    # Eq.: Briegleb 1992  alpha and gamma for direct radiation

    rdif_a[:, lyr] = (
        (ue**2 - 1) * (1 / extins - extins) / ne
    )  # R BAR = layer reflectivity to DIFFUSE radiation
    # T BAR layer transmissivity to DIFFUSE radiation
    tdif_a[:, lyr] = 4 * ue / ne

    # evaluate rdir, tdir for direct beam
    trnlay[:, lyr] = np.maximum(
        np.full((model_config.nbr_wvl,), exp_min), np.exp(-ts / mu0n)
    )  # transmission from TOA to interface

    #  Eq. 50: Briegleb and Light 2007  alpha and gamma for direct radiation
    alp = (
        (0.75 * ws * mu0n) * (1 + gs * (1 - ws)) / (1 - lm**2 * mu0n**2 + epsilon)
    )  # alp = alpha(ws,mu0n,gs,lm)
    gam = (0.5 * ws) * (
        (1 + 3 * gs * mu0n**2 * (1 - ws)) / (1 - lm**2 * mu0n**2 + epsilon)
    )  # gam = gamma(ws,mu0n,gs,lm)

    # apg = alpha plus gamma
    # amg = alpha minus gamma
    apg = alp + gam
    amg = alp - gam

    rdir[:, lyr] = apg * rdif_a[:, lyr] + amg * (
        tdif_a[:, lyr] * trnlay[:, lyr] - 1
    )  # layer reflectivity to DIRECT radiation
    tdir[:, lyr] = (
        apg * tdif_a[:, lyr] + (amg * rdif_a[:, lyr] - apg + 1) * trnlay[:, lyr]
    )  # layer transmissivity to DIRECT radiation

    return rdir, tdir, ts, ws, gs, lm


def define_constants_arrays(tau, g, ssa, illumination, ice, model_config):
    """defines and instantiates constants required for a-d calculations.

    Defines and instantiates all variables required for calculating energy fluxes
    using the adding-doubling method.

    Args:
        tau: optical thickness of ice column in m/m
        g: asymmetry parameter for ice column
        ssa: single scattering albedo of ice column (dimensionless)
        ice: instance of Ice class
        illumination: instance of Illumination class
        model_config: instance of ModelConfig class

    Returns:
        tau0: initial optical thickness (m/m)
        g0: initial asymmetry parameter (dimensionless)
        ssa0: initial single scatterign albedo (dimensionless)
        epsilon: small number to avoid singularity
        exp_min: small number to avoid zero calcs
        nr: real part of refractive index
        mu0: cosine of direct beam zenith angle
        mu0n: adjusted cosine of direct beam zenith angle after refraction
        trnlay: transmission through layer
        rdif_a: reflectivity to diffuse irradiance at polarization angle == perpendicular
        rdif_b: reflectivity to diffuse irradiance at polarization angle == parallel
        tdif_a: transmissivity to diffuse irradiance at polarization angle == perpendicular
        tdif_b: transmissivity to diffuse irradiance at polarization angle == parallel
        rdir: reflectivity to direct beam
        tdir: transmissivity to direct beam
        lyrfrsnl: index of uppermost fresnel reflecting layer in ice column
        trnlay:
        rdif_a:
        rdif_b:
        tdif_a:
        tdif_b:
        rdir:
        tdir:
        rdndif:
        trntdr:
        trndir:
        trndif:
        fdirup:
        fdifup:
        fdirdn:
        fdifdn:
        dfdir:
        dfdif:
        F_up:
        F_dwn:
        F_abs:
        F_abs_vis:
        F_abs_nir:
        rupdif:
        rupdir:

    """

    tau0 = tau.T  # read and transpose tau
    g0 = g.T  # read and transpose g
    ssa0 = ssa.T  # read and transpose ssa
    epsilon = 1e-5  # to deal with singularity
    exp_min = 1e-5  # exp(-500)  # min value > 0 to avoid error
    nr = np.zeros(shape=480)
    mu0 = illumination.mu_not * np.ones(480)  # cos beam angle = incident beam

    # ice-adjusted real refractive index
    temp1 = (
        ice.ref_idx_re**2
        - ice.ref_idx_im**2
        + np.sin(np.arccos(illumination.mu_not)) ** 2
    )
    temp2 = (
        ice.ref_idx_re**2
        - ice.ref_idx_im**2
        - np.sin(np.arccos(illumination.mu_not)) ** 2
    )
    nr = (np.sqrt(2) / 2) * (
        temp1 + (temp2**2 + 4 * ice.ref_idx_re**2 * ice.ref_idx_im**2) ** 0.5
    ) ** 0.5

    # . Eq. 20: Briegleb and Light 2007: adjusts beam angle
    # (i.e. this is Snell's Law for refraction at interface between media)
    # mu0n = -1 represents light travelling vertically upwards and mu0n = +1
    # represents light travellign vertically downwards
    # mu0n = np.sqrt(1-((1-mu0**2)/(ref_indx*ref_indx)))  (original,
    # before update for diffuse Fresnel reflection)
    # this version accounts for diffuse fresnel reflection:
    mu0n = np.cos(np.arcsin(np.sin(np.arccos(mu0)) / nr))

    # solar beam transm for layer (direct beam only)
    trnlay = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer reflectivity to diffuse radiation from above
    rdif_a = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer reflectivity to diffuse radiation from below
    rdif_b = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer transmission to diffuse radiation from above
    tdif_a = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer transmission to diffuse radiation from below
    tdif_b = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer reflectivity to direct radiation (solar beam + diffuse)
    rdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    # layer transmission to direct radiation (solar beam + diffuse)
    tdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])

    rdndif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trntdr = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trndif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trndir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trntdr[:, 0] = 1
    trndif[:, 0] = 1
    rdndif[:, 0] = 0
    trndir[:, 0] = 1

    fdirup = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdifup = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdirdn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdifdn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    dfdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    dfdif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    F_up = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    F_dwn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    F_abs = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr])
    F_abs_vis = np.zeros(shape=[ice.nbr_lyr])
    F_abs_nir = np.zeros(shape=[ice.nbr_lyr])

    # reflectivity to diffuse radiation
    rupdif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rupdif[:, ice.nbr_lyr] = ice.sfc
    # reflectivity to direct radiation
    rupdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rupdir[:, ice.nbr_lyr] = ice.sfc

    # if there are non zeros in layer type, grab the index of the
    # first fresnel layer and load in the precalculated diffuse fresnel
    # reflection
    # (precalculated as large no. of gaussian points required for convergence)
    if np.sum(np.array(ice.layer_type) == 1) > 0:
        lyrfrsnl = ice.layer_type.index(1)

    else:
        lyrfrsnl = 9999999

    return (
        tau0,
        g0,
        ssa0,
        epsilon,
        exp_min,
        nr,
        mu0,
        mu0n,
        trnlay,
        rdif_a,
        rdif_b,
        tdif_a,
        tdif_b,
        rdir,
        tdir,
        lyrfrsnl,
        trnlay,
        rdif_a,
        rdif_b,
        tdif_a,
        tdif_b,
        rdir,
        tdir,
        rdndif,
        trntdr,
        trndir,
        trndif,
        fdirup,
        fdifup,
        fdirdn,
        fdifdn,
        dfdir,
        dfdif,
        F_up,
        F_dwn,
        F_abs,
        F_abs_vis,
        F_abs_nir,
        rupdif,
        rupdir,
    )


def calc_reflection_transmission_from_top(
    lyr,
    trnlay,
    rdif_a,
    rdir,
    tdif_a,
    rdif_b,
    tdir,
    tdif_b,
    model_config,
    ice,
    trndir,
    rdndif,
    trntdr,
    trndif,
):
    """Calculates the reflection and transmission of energy at top surfaces.

    Calculate the solar beam transmission, total transmission, and
    reflectivity for diffuse radiation from below at interface lyr,
    the top of the current layer lyr.

    Args:
        lyr: integer representing the index of the current layer (0 at top)
        trnlay: transmissivity of current layer
        rdif_a: reflectivity to diffuse irradiance at polarization state == perpendicular
        rdir: reflectivity to direct beam
        tdif_a: transmissivity to diffuse irradiance at polarization state == perpendicular
        rdif_b: reflectivity to diffuse irradiance at polarization state == parallel
        tdir: transmissivity to direct beam
        tdif_b: transmissivity to diffuse irradiance at polarization state == parallel
        model_config: instance of ModelConfig class
        ice: instance of Ice class
        trndir:
        rdndif:
        trntdr:
        trndif:

    Returns:
        trndir: transmission of direct beam
        trntdr: total transmission of direct beam for all layers above current layer
        rdndif: downwards diffuse reflectance
        trndif: diffuse transmission

    """

    # Eq. 51  Briegleb and Light 2007
    trndir[:, lyr + 1] = (
        trndir[:, lyr] * trnlay[:, lyr]
    )  # solar beam transmission from top

    # interface multiple scattering for lyr-1
    refkm1 = 1 / (1 - rdndif[:, lyr] * rdif_a[:, lyr])

    # direct tran times layer direct ref
    tdrrdir = trndir[:, lyr] * rdir[:, lyr]

    # total down diffuse = tot tran - direct tran
    tdndif = trntdr[:, lyr] - trndir[:, lyr]

    # total transmission to direct beam for layers above
    trntdr[:, lyr + 1] = (
        trndir[:, lyr] * tdir[:, lyr]
        + (tdndif + tdrrdir * rdndif[:, lyr]) * refkm1 * tdif_a[:, lyr]
    )

    # Eq. B4  Briegleb and Light 2007
    # reflectivity to diffuse radiation for layers above
    rdndif[:, lyr + 1] = rdif_b[:, lyr] + (
        tdif_b[:, lyr] * rdndif[:, lyr] * refkm1 * tdif_a[:, lyr]
    )
    # diffuse transmission to diffuse beam for layers above
    trndif[:, lyr + 1] = trndif[:, lyr] * refkm1 * tdif_a[:, lyr]

    return trndir, trntdr, rdndif, trndif


def apply_gaussian_integral(
    model_config, exp_min, ts, ws, gs, epsilon, lm, lyr, rdif_a, tdif_a
):
    """Applies gaussian integral to integrate over angles.

    Uses gaussien integration to integrate fluxes hemispherically from
    N of reference angles where N = len(gauspt) (default is 8).

    Args:
        model_config: instance of ModelConfig class
        exp_min: small number for avoiding div/0 error
        ts: delta-scaled extinction optical depth for lyr
        ws: delta-scaled single scattering albedo for lyr
        gs: delta-scaled asymmetry parameter for lyr
        epsilon: small number to avoid singularity
        lm: lamda for use in delta scaling
        lyr: integer representing index of current layer (0==top)
        rdif_a: rdif_a: reflectance to diffuse energy w polarization state == perpendicular
        tdif_a: rdif_a: transmittance to diffuse energy w polarization state == perpendicular

    Returns:
        smt: accumulator for tdif gaussian integration
        smr: accumulator for rdif gaussian integration
        swt: sum of gaussian weights

    """
    # gaussian angles (radians)
    gauspt = [
        0.9894009,
        0.9445750,
        0.8656312,
        0.7554044,
        0.6178762,
        0.4580168,
        0.2816036,
        0.0950125,
    ]
    # gaussian weights
    gauswt = [
        0.0271525,
        0.0622535,
        0.0951585,
        0.1246290,
        0.1495960,
        0.1691565,
        0.1826034,
        0.1894506,
    ]
    swt = 0
    smr = 0
    smt = 0
    R1 = rdif_a[:, lyr]  # use R1 as temporary var
    T1 = tdif_a[:, lyr]  # use T1 as temporary var

    for ng in np.arange(0, len(gauspt), 1):
        mu = gauspt[ng]  # solar zenith angles
        gwt = gauswt[ng]  # gaussian weight
        swt = swt + mu * gwt  # sum of weights
        trn = np.maximum(
            np.full((model_config.nbr_wvl,), exp_min), np.exp(-ts / mu)
        )  # transmission

        alp = (
            (0.75 * ws * mu) * (1 + gs * (1 - ws)) / (1 - lm**2 * mu**2 + epsilon)
        )  # alp = alpha(ws,mu0n,gs,lm)
        gam = (
            (0.5 * ws)
            * (1 + 3 * gs * mu**2 * (1 - ws))
            / (1 - lm**2 * mu**2 + epsilon)
        )  # gam = gamma(ws,mu0n,gs,lm)

        apg = alp + gam
        amg = alp - gam
        rdr = apg * R1 + amg * T1 * trn - amg
        tdr = apg * T1 + amg * R1 * trn - apg * trn + trn
        smr = smr + mu * rdr * gwt  # accumulator for rdif gaussian integration
        smt = smt + mu * tdr * gwt  # accumulator for tdif gaussian integration

    return smt, smr, swt


def update_transmittivity_reflectivity(
    swt, smr, smt, lyr, rdif_a, tdif_a, rdif_b, tdif_b
):
    """updates transmissivity and reflectivity values after iterations.

    Args:
        swt: sum of gaussian weights (for integrating over angle)
        smr: accumulator for rdif gaussian integration
        smt: accumulator for tdif gaussian integration
        lyr: integer representign index of current layer (0 == top)
        rdif_a: reflectance to diffuse energy w polarization state == perpendicular
        rdif_b: reflectance to diffuse energy w polarization state == parallel
        tdif_a: transmittance to diffuse energy w polarization state == perpendicular
        tdif_b: transmittance to diffuse energy w polarization state == parallel

    Returns:
        rdif_a: updated reflectance to diffuse energy w polarization state == perpendicular
        rdif_b: updated reflectance to diffuse energy w polarization state == parallel
        tdif_a: updated transmittance to diffuse energy w polarization state == perpendicular
        tdif_b: updated transmittance to diffuse energy w polarization state == parallel
    """
    rdif_a[:, lyr] = smr / swt
    tdif_a[:, lyr] = smt / swt

    # homogeneous layer
    rdif_b[:, lyr] = rdif_a[:, lyr]
    tdif_b[:, lyr] = tdif_a[:, lyr]

    return rdif_a, tdif_a, rdif_b, tdif_b


def calc_correction_fresnel_layer(
    model_config,
    ice,
    illumination,
    mu0n,
    mu0,
    nr,
    rdif_a,
    rdif_b,
    tdif_a,
    tdif_b,
    trnlay,
    lyr,
    rdir,
    tdir,
):
    """Calculates correction for Fresnel reflection and total internal reflection.

    Corrects fluxes for Fresnel reflection in cases where total
    internal reflection does and does not occur (angle > critical_angle).
    In TIR case fluxes are precalculated because ~256 gaussian points required
    for convergence.

    Args:
        model_config: instance of ModelConfig class
        ice: instance of Ice class
        illumination: instance of Illumination class
        mu0n: incidence angle for direct beam adjusted for refraction
        mu0: incidence angle of direct beam at upper surface
        nr: real part of refractive index
        rdif_a: reflectance to diffuse energy w polarization state == perpendicular
        rdif_b: reflectance to diffuse energy w polarization state == parallel
        tdif_a: transmittance to diffuse energy w polarization state == perpendicular
        tdif_b: transmittance to diffuse energy w polarization state == parallel
        trnlay: transmission of layer == lyr
        lyr: current layer (0 ==top)
        rdir: reflectance to direct beam
        tdir: transmission of direct beam


    Returns:
        rdif_a: updated reflectance to diffuse energy w polarization state == perpendicular
        rdif_b: updated reflectance to diffuse energy w polarization state == parallel
        tdif_a: updated transmittance to diffuse energy w polarization state == perpendicular
        tdif_b: updated transmittance to diffuse energy w polarization state == parallel
        trnlay: updated transmission of layer == lyr
        rdir: updated reflectance to direct beam
        tdir: updated transmission of direct beam
    """

    ref_indx = ice.ref_idx_re + 1j * ice.ref_idx_im
    critical_angle = np.arcsin(ref_indx)

    for wl in np.arange(0, model_config.nbr_wvl, 1):
        if np.arccos(illumination.mu_not) < critical_angle[wl]:
            # in this case, no total internal reflection

            # compute fresnel reflection and transmission amplitudes
            # for two polarizations: 1=perpendicular and 2=parallel to
            # the plane containing incident, reflected and refracted rays.

            # Eq. 22  Briegleb & Light 2007
            # Inputs to equation 21 (i.e. Fresnel formulae for R and T)
            R1 = (mu0[wl] - nr[wl] * mu0n[wl]) / (
                mu0[wl] + nr[wl] * mu0n[wl]
            )  # reflection amplitude factor for perpendicular polarization
            R2 = (nr[wl] * mu0[wl] - mu0n[wl]) / (
                nr[wl] * mu0[wl] + mu0n[wl]
            )  # reflection amplitude factor for parallel polarization
            T1 = (
                2 * mu0[wl] / (mu0[wl] + nr[wl] * mu0n[wl])
            )  # transmission amplitude factor for perpendicular polarization
            T2 = (
                2 * mu0[wl] / (nr[wl] * mu0[wl] + mu0n[wl])
            )  # transmission amplitude factor for parallel polarization

            # unpolarized light for direct beam
            # Eq. 21  Brigleb and light 2007
            Rf_dir_a = 0.5 * (R1**2 + R2**2)
            Tf_dir_a = 0.5 * (T1**2 + T2**2) * nr[wl] * mu0n[wl] / mu0[wl]

        else:  # in this case, total internal reflection occurs
            Tf_dir_a = 0
            Rf_dir_a = 1

        # precalculated diffuse reflectivities and transmissivities
        # for incident radiation above and below fresnel layer, using
        # the direct albedos and accounting for complete internal
        # reflection from below. Precalculated because high order
        # number of gaussian points (~256) is required for convergence:

        # Eq. 25  Brigleb and light 2007
        # diffuse reflection of flux arriving from above

        # reflection from diffuse unpolarized radiation
        Rf_dif_a = ice.fl_r_dif_a[wl]
        Tf_dif_a = 1 - Rf_dif_a  # transmission from diffuse unpolarized radiation

        # diffuse reflection of flux arriving from below
        Rf_dif_b = ice.fl_r_dif_b[wl]
        Tf_dif_b = 1 - Rf_dif_b

        # -----------------------------------------------------------------------
        # the lyr = lyrfrsnl layer properties are updated to combine
        # the fresnel (refractive) layer, always taken to be above
        # the present layer lyr (i.e. be the top interface):

        # denom interface scattering
        rintfc = 1 / (1 - Rf_dif_b * rdif_a[wl, lyr])

        # layer transmissivity to DIRECT radiation
        # Eq. B7  Briegleb & Light 2007
        tdir[wl, lyr] = (
            Tf_dir_a * tdir[wl, lyr]
            + Tf_dir_a * rdir[wl, lyr] * Rf_dif_b * rintfc * tdif_a[wl, lyr]
        )

        # layer reflectivity to DIRECT radiation
        # Eq. B7  Briegleb & Light 2007
        rdir[wl, lyr] = Rf_dir_a + Tf_dir_a * rdir[wl, lyr] * rintfc * Tf_dif_b

        # R BAR = layer reflectivity to DIFFUSE radiation (above)
        # Eq. B9  Briegleb & Light 2007
        rdif_a[wl, lyr] = Rf_dif_a + Tf_dif_a * rdif_a[wl, lyr] * rintfc * Tf_dif_b

        # R BAR = layer reflectivity to DIFFUSE radiation (below)
        # Eq. B10  Briegleb & Light 2007
        rdif_b[wl, lyr] = (
            rdif_b[wl, lyr] + tdif_b[wl, lyr] * Rf_dif_b * rintfc * tdif_a[wl, lyr]
        )

        # T BAR layer transmissivity to DIFFUSE radiation (above),
        # Eq. B9  Briegleb & Light 2007
        tdif_a[wl, lyr] = tdif_a[wl, lyr] * rintfc * Tf_dif_a

        # Eq. B10  Briegleb & Light 2007
        tdif_b[wl, lyr] = tdif_b[wl, lyr] * rintfc * Tf_dif_b

        # update trnlay to include fresnel transmission
        trnlay[wl, lyr] = Tf_dir_a * trnlay[wl, lyr]

    return rdif_a, rdif_b, tdif_a, tdif_b, trnlay, rdir, tdir


def calc_reflection_below(
    ice,
    model_config,
    rdif_a,
    rdif_b,
    tdif_a,
    tdif_b,
    trnlay,
    rdir,
    tdir,
    rupdif,
    rupdir,
):
    """Calculates dir/diff reflectyivity for layers below surface.

    Compute reflectivity to direct (rupdir) and diffuse (rupdif) radiation
    for layers below by adding succesive layers starting from the
    underlying ice and working upwards.

    Args:
        model_config: instance of ModelConfig class
        ice: instance of Ice class
        rdif_a: reflectance to diffuse energy w polarization state == perpendicular
        rdif_b: reflectance to diffuse energy w polarization state == parallel
        tdif_a: transmittance to diffuse energy w polarization state == perpendicular
        tdif_b: transmittance to diffuse energy w polarization state == parallel
        trnlay: transmission of layer == lyr
        rdir: reflectance to direct beam
        tdir: transmission of direct beam
        rupdif: upwards flux direct
        rupdir: upwards flux diffuse

    Returns:
        rupdir: upwards flux direct
        rupdif: upwards flux diffuse

    """

    for lyr in np.arange(
        ice.nbr_lyr - 1, -1, -1
    ):  # starts at the bottom and works its way up to the top layer
        # Eq. B5  Briegleb and Light 2007
        # interface scattering
        refkp1 = 1 / (1 - rdif_b[:, lyr] * rupdif[:, lyr + 1])

        # dir from top layer plus exp tran ref from lower layer, interface
        # scattered and tran thru top layer from below, plus diff tran ref
        # from lower layer with interface scattering tran thru top from below
        rupdir[:, lyr] = (
            rdir[:, lyr]
            + (
                trnlay[:, lyr] * rupdir[:, lyr + 1]
                + (tdir[:, lyr] - trnlay[:, lyr]) * rupdif[:, lyr + 1]
            )
            * refkp1
            * tdif_b[:, lyr]
        )

        # dif from top layer from above, plus dif tran upwards reflected and
        # interface scattered which tran top from below
        rupdif[:, lyr] = (
            rdif_a[:, lyr]
            + tdif_a[:, lyr] * rupdif[:, lyr + 1] * refkp1 * tdif_b[:, lyr]
        )
    return rupdif, rupdir


def trans_refl_at_interfaces(
    model_config,
    ice,
    rupdif,
    rupdir,
    rdndif,
    trndir,
    trndif,
    trntdr,
    fdirup,
    fdirdn,
    fdifup,
    fdifdn,
    dfdir,
    dfdif,
):
    """Calculates transmission and reflection at layer interfaces.

    Args:
        model_config: instance of ModelConfig class
        ice: instance of Ice class
        rupdif: total diffuse radiation reflected upwards
        rupdir: total direct radiation reflected upwards
        rdndif: downwards reflection of diffuse radiation
        trndir: transmission of direct radiation
        trndif: transmission of diffuse radiation
        trntdr: total transmission
        fdirup:
        fdirdn:
        fdifup:
        fdifdn:
        dfdir:
        dfdif:

    Returns:
        fdirup: upwards flux of direct radiation
        fdifup: upwards flux of diffuse radiation
        fdirdn: downwards flux of direct radiation
        fdifdn: downwards flux of diffuse radiation

    """

    puny = 1e-10  # not sure how should we define this

    for lyr in np.arange(0, ice.nbr_lyr + 1, 1):
        # Eq. 52  Briegleb and Light 2007
        # interface scattering
        refk = 1 / (1 - rdndif[:, lyr] * rupdif[:, lyr])

        # dir tran ref from below times interface scattering, plus diff
        # tran and ref from below times interface scattering
        fdirup[:, lyr] = (
            trndir[:, lyr] * rupdir[:, lyr]
            + (trntdr[:, lyr] - trndir[:, lyr]) * rupdif[:, lyr]
        ) * refk

        # dir tran plus total diff trans times interface scattering plus
        # dir tran with up dir ref and down dif ref times interface scattering
        fdirdn[:, lyr] = (
            trndir[:, lyr]
            + (
                trntdr[:, lyr]
                - trndir[:, lyr]
                + trndir[:, lyr] * rupdir[:, lyr] * rdndif[:, lyr]
            )
            * refk
        )

        # diffuse tran ref from below times interface scattering
        fdifup[:, lyr] = trndif[:, lyr] * rupdif[:, lyr] * refk

        # diffuse tran times interface scattering
        fdifdn[:, lyr] = trndif[:, lyr] * refk

        # dfdir = fdirdn - fdirup
        dfdir[:, lyr] = (
            trndir[:, lyr]
            + (trntdr[:, lyr] - trndir[:, lyr]) * (1 - rupdif[:, lyr]) * refk
            - trndir[:, lyr] * rupdir[:, lyr] * (1 - rdndif[:, lyr]) * refk
        )

        if np.max(dfdir[:, lyr]) < puny:
            dfdir[:, lyr] = np.zeros(
                (model_config.nbr_wvl,), dtype=int
            )  # echmod necessary?
            # dfdif = fdifdn - fdifup

        dfdif[:, lyr] = trndif[:, lyr] * (1 - rupdif[:, lyr]) * refk

        if np.max(dfdif[:, lyr]) < puny:
            dfdif[:, lyr] = np.zeros(
                (model_config.nbr_wvl,), dtype=int
            )  # !echmod necessary?

    return fdirup, fdifup, fdirdn, fdifdn


def calculate_fluxes(
    model_config,
    ice,
    illumination,
    fdirup,
    fdifup,
    fdirdn,
    fdifdn,
    F_up,
    F_dwn,
    F_abs,
    F_abs_vis,
    F_abs_nir,
):
    """Calculates total fluxes in each layer and for entire column.

    Args:
        model_config: instance of ModelConfig class
        ice: instance of Ice class
        illumination: instance of Illumination class
        fdirup: upwards flux of direct radiation
        fdifup: upwards flux od diffuse radiation
        fdirdn: downwards flux of direct radiation
        fdifdn: downwards flux of diffuse radiation
        F_up:
        F_dwn:
        F_abs:
        F_abs_vis:
        F_abs_nir:

    Returns:
        albedo: ratio of upwards fluxes to incoming irradiance
        F_abs: absorbed flux in each layer
        F_btm_net: net fluxes at bottom surface
        F_top_pls: upwards flux from upper surface

    """

    for n in np.arange(0, ice.nbr_lyr + 1, 1):
        F_up[:, n] = (
            fdirup[:, n] * (illumination.Fs * illumination.mu_not * np.pi)
            + fdifup[:, n] * illumination.Fd
        )
        F_dwn[:, n] = (
            fdirdn[:, n] * (illumination.Fs * illumination.mu_not * np.pi)
            + fdifdn[:, n] * illumination.Fd
        )

    F_net = F_up - F_dwn

    # import matplotlib.pyplot as plt
    # plt.plot(Inputs.Fs)
    # plt.show()
    # plt.figure(1)
    # plt.plot(Inputs.Fs)
    # plt.plot(F_dwn[:,0].T)
    # plt.plot(F_dwn[:,1].T)
    # plt.plot(F_dwn[:,2].T)
    # plt.show()
    # plt.show()

    # Absorbed flux in each layer
    F_abs[:, :] = F_net[:, 1:] - F_net[:, :-1]

    # Upward flux at upper model boundary
    F_top_pls = F_up[:, 0]

    # Net flux at lower model boundary = bulk transmission through entire
    # media = absorbed radiation by underlying surface:
    F_btm_net = -F_net[:, ice.nbr_lyr]

    for i in np.arange(0, ice.nbr_lyr, 1):  # [0,1,2,3,4]
        F_abs_vis[i] = sum(F_abs[0 : model_config.vis_max_idx, i])
        F_abs_nir[i] = sum(
            F_abs[model_config.vis_max_idx : model_config.nir_max_idx, i]
        )

    albedo = F_up[:, 0] / F_dwn[:, 0]

    return albedo, F_abs, F_btm_net, F_top_pls


def conservation_of_energy_check(illumination, F_abs, F_btm_net, F_top_pls):
    """Checks there is no conservation of energy violation.

    Args:
        illumination: instance of Illumination class
        F_abs: absorbed flux in each layer
        F_btm_net: net flux at bottom surface
        F_top_pls: upwards flux from upper boundary

    Returns:
        None

    Raises:
        ValueError is conservation of energy error is detected

    """
    # Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
    energy_sum = (
        (illumination.mu_not * np.pi * illumination.Fs)
        + illumination.Fd
        - (np.sum(F_abs, axis=1) + F_btm_net + F_top_pls)
    )

    energy_conservation_error = sum(abs(energy_sum))

    if energy_conservation_error > 1e-10:
        raise ValueError(f"energy conservation error: {energy_conservation_error}")
    else:
        pass


def get_outputs(illumination, albedo, model_config, L_snw, F_abs, F_btm_net):
    """Assimilates useful data into instance of Outputs class.

    Args:
        illumination: instance of Illumination class
        albedo: ratio of upwwards fluxes and irradiance
        model_config: instance of ModelConfig class
        L_snw: mass of ice in each layer
        F_abs: absorbed flux in each layer
        F_btm_net: net flux at bottom surface

    Returns:
        outputs: instance of Outputs class

    """
    outputs = Outputs()

    # Radiative heating rate:
    F_abs_slr = np.sum(F_abs, axis=0)
    # [K/s] 2117 = specific heat ice (J kg-1 K-1)
    heat_rt = F_abs_slr / (L_snw * 2117)
    outputs.heat_rt = heat_rt * 3600  # [K/hr]

    # Spectral albedo
    outputs.albedo = albedo

    # Spectrally-integrated solar, visible, and NIR albedos:
    outputs.BBA = np.sum(illumination.flx_slr * albedo) / np.sum(illumination.flx_slr)
    outputs.BBAVIS = np.sum(
        illumination.flx_slr[0 : model_config.vis_max_idx]
        * albedo[0 : model_config.vis_max_idx]
    ) / np.sum(illumination.flx_slr[0 : model_config.vis_max_idx])
    outputs.BBANIR = np.sum(
        illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx]
        * albedo[model_config.vis_max_idx : model_config.nir_max_idx]
    ) / np.sum(
        illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx]
    )

    # Total incident insolation( Wm - 2)
    outputs.total_insolation = np.sum(
        (illumination.mu_not * np.pi * illumination.Fs) + illumination.Fd
    )

    # Spectrally-integrated absorption by underlying surface:
    outputs.abs_slr_btm = np.sum(F_btm_net, axis=0)
    outputs.abs_vis_btm = np.sum(F_btm_net[0 : model_config.vis_max_idx], axis=0)
    outputs.abs_nir_btm = np.sum(
        F_btm_net[model_config.vis_max_idx : model_config.nir_max_idx + 1], axis=0
    )

    # Spectrally-integrated VIS and NIR total snowpack absorption:
    outputs.abs_vis_tot = np.sum(
        illumination.flx_slr[0 : model_config.vis_max_idx]
        * (1 - albedo[0 : model_config.vis_max_idx])
    )
    outputs.abs_nir_tot = np.sum(
        illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx]
        * (1 - albedo[model_config.vis_max_idx : model_config.nir_max_idx])
    )
    # Spectrally-integrated absorption by entire snow/ice column
    outputs.abs_slr_tot = np.sum(F_abs_slr)

    # Spectrally-integrated absorption by each layer
    outputs.absorbed_flux_per_layer = F_abs_slr

    return outputs


def apply_smoothing_function(albedo, model_config):
    """Applies Savitsky-Golay smoothing function to albedo, if toggled.

    Args:
        albedo: array of albedo values, likely passed as outputs.albedo
        model_config: instance of ModelConfig

    Returns:
        albedo: updated array of albedo values

    """

    yhat = savgol_filter(albedo, model_config.window_size, model_config.poly_order)
    albedo = yhat

    return albedo


if __name__ == "__main__":
    pass
