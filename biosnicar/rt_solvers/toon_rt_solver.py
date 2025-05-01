#!/usr/bin/python
"""Radiative transfer solver using the Toon matrix method.

Here the ice/impurity optical properties and illumination conditions
are used to calculate energy fluxes between the ice, atmosphere and
underlying substrate.

Typically this function would be called from snicar_driver() because
it takes as inputs intermediates that are calculated elsewhere.
Specifically, the functions setup_snicar(), get_layer_OPs() and
mix_in_impurities() are called to generate tau, ssa, g and L_snw,
which are then passed as inputs to toon_solver().

The Toon solver was originally coded in FORTRAN and then Matlab by Mark Flanner.
His Matlab scripts were the jumping off point for this script and are still
used to benchmark this code against.

This solver can only treat ice as a "granular" material with a bulk medium of air
with discrete ice grains. No Fresnel reflection is included in this solver. There
are also well-known instabilities that arise at certain solar zenith angles (see
Toon et al. 1989). Our model does not accept solar zenith angles that fall outside
of the valid range for the Toon solver. Swapping to the adding-doubling solver is
a route around this problem - there the solar zenith angle can range from 1-89
degrees. The main benefit of the Toon solver is speed, and continuity from legacy
radiative transfer models.

"""
import numpy as np
from scipy.signal import savgol_filter

from biosnicar.classes.outputs import Outputs


def toon_solver(tau, ssa, g, L_snw, ice, illumination, model_config, rt_config):
    """Driver func for toon radiative transfer solver.

    Makes function calls relating to radiative transfer solver in sequence and returns outputs.

    Args:
        tau: optical thickness
        ssa: single scattering albedo
        g: asymmetry parameter
        L_snw: mass of ice in each layer
        ice: instance of Ice class
        illumination: instance of Illumination class
        model_config: instance of ModelConfig class
        rt_config: instance of RTConfig class

    Returns:
        outputs: instance of Outputs class

    """
    validate_inputs_toon(ice, illumination)
    # ----------------------------------------------------------------------------------
    # PERFORM DELTA TRANSFORMATION IF REQUIRED
    # ----------------------------------------------------------------------------------
    # The star represents the delta transformed quantity
    # if no delta transformation is applied, the starred quantity
    # is equal to the unstarred quantity

    g_star, ssa_star, tau_star = delta_transformation(rt_config, g, ssa, tau)

    # CALCULATE TOTAL OPTICAL DEPTH OF ENTIRE COLUMN
    # i.e. tau_clm = total optical depth from upper boundary
    # to upper boundary of layer n. This is therefore a cumulative
    # quantity - subsequently lower layers contain the sum of the
    # # optical depth of all overlying layers

    tau_clm = calculate_optical_depth_of_column(ice, model_config, tau_star)

    # SET BOUNDARY CONDITION: BOTTOM BOUNDARY
    # calculate radiation reflected skywards by underlying surface
    # (i.e. lower model boundary)
    # remainder is lost

    s_sfc = boundary_condition(ice, illumination, tau_clm, tau_star)

    # ----------------------------------------------------------------------------------
    # Apply Two-Stream Approximation (Toon et al, table 1)
    # ----------------------------------------------------------------------------------
    gamma1, gamma2, gamma3, gamma4, mu_one = two_stream_approximation(
        rt_config, ssa_star, g_star, illumination
    )

    # Toon et al equation 21 and 22
    # Note that the values of lam and GAMMA depend upon gamma1 and gamma2, which
    # vary depending upon the two-stream approximation used
    # variable "lambda" renamed "lam" to avoid confusion with lambda function
    lam, GAMMA, e1, e2, e3, e4 = calculate_matrix_coefficients(gamma1, gamma2, tau_star)

    # ----------------------------------------------------------------------------------
    # Calculate C-functions
    # ----------------------------------------------------------------------------------

    # C is the direct beam flux calculated at the top and bottom of each layer, i,
    # see Toon equations 23 and 24
    C_pls_btm, C_mns_btm, C_pls_top, C_mns_top = c_functions(
        ice,
        model_config,
        illumination,
        ssa_star,
        tau_star,
        tau_clm,
        lam,
        gamma1,
        gamma2,
        gamma3,
        gamma4,
    )

    # ----------------------------------------------------------------------------------
    # Initialize tridiagonal matrix solution
    # ----------------------------------------------------------------------------------

    # expanding the number of layers to 2*ice.nbr_lyr so that fluxes at upper and lower
    # layer boundaries can be resolved. This section was confusing to code -
    # for each layer
    # index (n) a second pair of indices (2 x i) are required. Different solutions are
    # applied depending upon whether i is even or odd. To translate the indexing for
    # this from FORTRAN/MATLAB into Python, it was necessary to assert n = (i/2)-1
    # for even layers and n = floor(i/2) for odd layers, with specific rules for the
    # boundaries i = 0 and i = ice.nbr_lyrs-1 (i.e. top surface and bottom surface).
    Y = matrix_solver(
        ice,
        illumination,
        model_config,
        s_sfc,
        C_pls_btm,
        C_mns_btm,
        C_pls_top,
        C_mns_top,
        e1,
        e2,
        e3,
        e4,
    )

    F_net, F_top_pls, F_btm_net, F_top_net, intensity = layer_fluxes(
        ice,
        illumination,
        model_config,
        tau_clm,
        tau_star,
        lam,
        Y,
        GAMMA,
        C_pls_btm,
        C_mns_btm,
        C_pls_top,
        C_mns_top,
        e1,
        e2,
        e3,
        e4,
        mu_one,
    )

    F_abs = absorbed_fluxes(ice, model_config, F_net, F_top_net)
    # Energy conservation check:
    # % Incident direct + diffuse radiation equals(absorbed + transmitted +
    # bulk_reflected)
    # conservation_of_energy_check(illumination, F_abs, F_btm_net, F_top_pls)
    outputs = get_outputs(
        model_config, ice, illumination, F_top_pls, F_top_net, F_btm_net, F_abs, L_snw
    )

    if model_config.smooth:
        outputs.albedo = apply_smoothing_function(outputs.albedo, model_config)

    return outputs


def validate_inputs_toon(ice, illumination):

    """Checks input parameters are valid for Toon solver.

    Checks for known invalid data configurations that break the Toon solver.
    If invalid config is detected a ValueError is raised, halting execution.

    Args:
        ice: class representing the ice column and containing related physical constants
        illumination: class representing incoming irradiance containing related physical constants

    Returns:
        None

    Raises:
        ValueError: raised with descriptive error message if invalid input detected
    """
    if np.sum(ice.layer_type) > 0:
        raise ValueError("There are solid ice layers in this model - use AD solver")

    elif (illumination.solzen < 50) or (illumination.solzen > 89):
        raise ValueError("Zenith angle out of valid range for Toon solver")

    elif np.sum(ice.cdom) > 0:
        raise ValueError("cdom is only available for solid ice layers")

    return


def delta_transformation(rt_config, g, ssa, tau):
    """Applied Delta transformation.

    Args:
        rt_config: instance of RTConfig class
        g: asymmetry parameter
        ssa: single scatterign albedo
        tau: optical thickness

    Returns:
        g_star: delta scaled asymmetry parameter
        ssa_star: delta scaled singloe scattering albedo
        tau_star: delta scaled optical thickness

    """

    if rt_config.delta:

        g_star = g / (1 + g)
        ssa_star = ((1 - (g**2)) * ssa) / (1 - (ssa * (g**2)))
        tau_star = (1 - (ssa * (g**2))) * tau

    else:

        g_star = g
        ssa_star = ssa
        tau_star = tau

    return g_star, ssa_star, tau_star


def calculate_optical_depth_of_column(ice, model_config, tau_star):
    """Calculates colum,n optical thickness.

    Args:
        ice: instance of Ice class
        model_config: instance of ModelConfig class
        tau_star: delta scaled optical thickness

    Returns:
        tau_clm: optical thickness of column

    """
    tau_clm = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])

    for i in np.arange(1, ice.nbr_lyr, 1):
        # start loop from 2nd layer, i.e. index = 1
        tau_clm[i, :] = tau_clm[i - 1, :] + tau_star[i - 1, :]

    return tau_clm


def boundary_condition(ice, illumination, tau_clm, tau_star):
    """Calculates reflectance from underlying surface.

    Args:
        ice: instance of Ice class
        illumination: instance of Illumination class
        tau_clm: optical thickness of column
        tau_star: delta_scaled optical thickness

    Returns:
        S_sfc: reflectance from underlying surface
    """
    s_sfc = (
        ice.sfc
        * illumination.mu_not
        * np.exp(
            -(tau_clm[ice.nbr_lyr - 1, :] + tau_star[ice.nbr_lyr - 1, :])
            / illumination.mu_not
        )
        * np.pi
        * illumination.Fs
    )
    return s_sfc


def two_stream_approximation(rt_config, ssa_star, g_star, illumination):
    """Applies two-stream approximation.

    Three 2-stream approximations are available: Eddington,
    Quadrature and hemispheric mean. The equations for each
    approximation are provided in Toon et al. (1989) Table 1.

    The hemispheric mean scheme is derived by assuming that the
    phase function is equal to 1  + g  in the forward scattering
    hemisphere and to 1  - g  in the backward scattering hemisphere.
    The asymmetry parameter is g. The hemispheric mean is only
    useful for infrared wavelengths.

    Args:
        rt_config: instance of RTConfig class
        ssa_star: delta scaled signle scattering albedo
        g_star: delta scaled asymmetry parameter

    Returns:
        gamma1: coefficient for matrix solution
        gamma2: coefficient for matrix solution
        gamma3: coefficient for matrix solution
        gamma4: coefficient for matrix solution
        mu_one: adjusted incidence angle

    """

    if rt_config.aprx_typ == 1:
        # apply Eddington approximation
        gamma1 = (7 - (ssa_star * (4 + (3 * g_star)))) / 4
        gamma2 = -(1 - (ssa_star * (4 - (3 * g_star)))) / 4
        gamma3 = (2 - (3 * g_star * illumination.mu_not)) / 4
        gamma4 = 1 - gamma3
        mu_one = 0.5

    elif rt_config.aprx_typ == 2:
        # apply quadrature approximation
        gamma1 = np.sqrt(3) * (2 - (ssa_star * (1 + g_star))) / 2
        gamma2 = ssa_star * np.sqrt(3) * (1 - g_star) / 2
        gamma3 = (1 - (np.sqrt(3) * g_star * illumination.mu_not)) / 2
        gamma4 = 1 - gamma3
        mu_one = 1 / np.sqrt(3)

    elif rt_config.aprx_typ == 3:
        # apply hemispheric mean approximation
        gamma1 = 2 - (ssa_star * (1 + g_star))
        gamma2 = ssa_star * (1 - g_star)
        gamma3 = (1 - (np.sqrt(3) * g_star * illumination.mu_not)) / 2
        gamma4 = 1 - gamma3
        mu_one = 0.5

    return gamma1, gamma2, gamma3, gamma4, mu_one


def calculate_matrix_coefficients(gamma1, gamma2, tau_star):
    """calculates coefficients required for matrix calculation.

    Args:
        gamma1: coefficient for matrix solution
        gamma2: coefficient for matrix solution
        tau_star: delta scaled optical thickness

    Returns:
        lam: coefficient for matrix solution
        GAMMA: coefficient for matrix solution
        e1: coefficient for matrix solution
        e2: coefficient for matrix solution
        e3: coefficient for matrix solution
        e4: coefficient for matrix solution
    """

    lam = np.sqrt(abs((gamma1**2) - (gamma2**2)))
    GAMMA = gamma2 / (gamma1 + lam)

    # calculate coefficients required for tridiagonal matrix calculation
    # (Toon et al Equation 44)
    e1 = 1 + (GAMMA * np.exp(-lam * tau_star))
    e2 = 1 - (GAMMA * np.exp(-lam * tau_star))
    e3 = GAMMA + np.exp(-lam * tau_star)
    e4 = GAMMA - np.exp(-lam * tau_star)

    return lam, GAMMA, e1, e2, e3, e4


def c_functions(
    ice,
    model_config,
    illumination,
    ssa_star,
    tau_star,
    tau_clm,
    lam,
    gamma1,
    gamma2,
    gamma3,
    gamma4,
):
    """Calculate fluxes through column.

    Args:
        ice: instance of Ice class
        model_config: instance of ModelConfig class
        illumination: instance of Illumination class
        ssa_star: delta scaled single scatterign albedo
        tau_star: delta scaled optical thickness
        tau_clm: column optical thickness
        lam: coefficient for matrix solution
        gamma1: coefficient for matrix solution
        gamma2: coefficient for matrix solution
        gamma3: coefficient for matrix solution
        gamma4: coefficient for matrix solution

    Returns:
        C_pls_btm: upwards flux from bottom of layer
        C_mns_btm: downwards flux from bottom of layer
        C_pls_top: upwards flux from top of layer
        C_mns_top: downwards flux from top of layer

    """
    C_pls_btm = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    C_mns_btm = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    C_pls_top = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    C_mns_top = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])

    for i in np.arange(0, ice.nbr_lyr, 1):

        if np.sum(illumination.Fs) > 0.0:

            C_pls_btm[i, :] = (
                ssa_star[i, :]
                * np.pi
                * illumination.Fs
                * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / illumination.mu_not)
                * (
                    ((gamma1[i, :] - (1 / illumination.mu_not)) * gamma3[i, :])
                    + (gamma4[i, :] * gamma2[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / (illumination.mu_not**2)))

            C_mns_btm[i, :] = (
                ssa_star[i, :]
                * np.pi
                * illumination.Fs
                * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / illumination.mu_not)
                * (
                    ((gamma1[i, :] + (1 / illumination.mu_not)) * gamma4[i, :])
                    + (gamma2[i, :] * gamma3[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / illumination.mu_not**2))

            C_pls_top[i, :] = (
                ssa_star[i, :]
                * np.pi
                * illumination.Fs
                * np.exp(-tau_clm[i, :] / illumination.mu_not)
                * (
                    (gamma1[i, :] - (1 / illumination.mu_not)) * gamma3[i, :]
                    + (gamma4[i, :] * gamma2[i, :])
                )
            ) / ((lam[1, :] ** 2) - (1 / illumination.mu_not**2))

            C_mns_top[i, :] = (
                ssa_star[i, :]
                * np.pi
                * illumination.Fs
                * np.exp(-tau_clm[i, :] / illumination.mu_not)
                * (
                    (gamma1[i, :] + (1 / illumination.mu_not)) * gamma4[i, :]
                    + (gamma2[i, :] * gamma3[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / illumination.mu_not**2))

        else:
            # no direct-beam flux:
            C_pls_btm[i, :] = 0
            C_mns_btm[i, :] = 0
            C_pls_top[i, :] = 0
            C_mns_top[i, :] = 0

    return C_pls_btm, C_mns_btm, C_pls_top, C_mns_top


def matrix_solver(
    ice,
    illumination,
    model_config,
    s_sfc,
    C_pls_btm,
    C_mns_btm,
    C_pls_top,
    C_mns_top,
    e1,
    e2,
    e3,
    e4,
):
    """Execute matrix calculation.

    Args:
        ice: instance of Ice class
        illumination: instance of Illumination class
        model_config: instance of ModelConfig class
        s_sfc: reflectance of underlying surface
        C_pls_btm: upwards flux from bottom of layer
        C_mns_btm: downwards flux from bottom of layer
        C_pls_top: upwards flux from top of layer
        C_mns_top: downwards flux from top of layer
        e1: coefficient for matrix solution
        e2: coefficient for matrix solution
        e3: coefficient for matrix solution
        e4: coefficient for matrix solution

    Returns:
        Y: intermediate representing t-0boundary fluxes
    """
    # Toon equations 41-43.
    # Boundary values for i=1 and i=2 * nbr_lyr, specifics for i=odd and i=even
    # Set up lists
    A = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])
    B = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])
    D = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])
    E = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])

    for i in np.arange(0, 2 * ice.nbr_lyr, 1):

        # TOP LAYER
        if i == 0:
            A[0, :] = 0.0
            B[0, :] = e1[0, :]
            D[0, :] = -e2[0, :]
            E[0, :] = illumination.Fd - C_mns_top[0, :]

        # BOTTOM LAYER
        elif i == 2 * ice.nbr_lyr - 1:
            A[i, :] = e1[ice.nbr_lyr - 1, :] - (ice.sfc * e3[ice.nbr_lyr - 1, :])
            B[i, :] = e2[ice.nbr_lyr - 1, :] - (ice.sfc * e4[ice.nbr_lyr - 1, :])
            D[i, :] = 0.0
            E[i, :] = (
                s_sfc[:]
                - C_pls_btm[ice.nbr_lyr - 1, :]
                + (ice.sfc * C_mns_btm[ice.nbr_lyr - 1, :])
            )

        # EVEN NUMBERED LAYERS
        elif i % 2 == 0:
            n = int(i / 2) - 1
            A[i, :] = (e2[n, :] * e3[n, :]) - (e4[n, :] * e1[n, :])
            B[i, :] = (e1[n, :] * e1[n + 1, :]) - (e3[n, :] * e3[n + 1, :])
            D[i, :] = (e3[n, :] * e4[n + 1, :]) - (e1[n, :] * e2[n + 1, :])
            E[i, :] = (e3[n, :] * (C_pls_top[n + 1, :] - C_pls_btm[n, :])) + (
                e1[n, :] * (C_mns_btm[n, :] - C_mns_top[n + 1, :])
            )

        # ODD NUMBERED LAYERS
        elif (i % 2 == 1) and (i < 2 * ice.nbr_lyr - 1):

            n = int(np.floor(i / 2))
            A[i, :] = (e2[n + 1, :] * e1[n, :]) - (e3[n, :] * e4[n + 1, :])
            B[i, :] = (e2[n, :] * e2[n + 1, :]) - (e4[n, :] * e4[n + 1, :])
            D[i, :] = (e1[n + 1, :] * e4[n + 1, :]) - (e2[n + 1, :] * e3[n + 1, :])
            E[i, :] = (e2[n + 1, :] * (C_pls_top[n + 1, :] - C_pls_btm[n, :])) + (
                e4[n + 1, :] * (C_mns_top[n + 1, :] - C_mns_btm[n, :])
            )

    # Now the actual tridiagonal matrix solving. Simply dividing A/B and E/B
    # throws an exception due to division by zero. Here we use numpy's nan_to_num
    # function to achieve the division where possible and replace nans with zeros.
    # We also set numpy to ignore the division error.

    # for bottom layer only
    # Toon et al Eq 45
    AS = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])
    DS = np.zeros([2 * ice.nbr_lyr, model_config.nbr_wvl])

    np.seterr(divide="ignore", invalid="ignore")
    AS[2 * ice.nbr_lyr - 1, :] = np.nan_to_num(
        A[2 * ice.nbr_lyr - 1, :] / B[2 * ice.nbr_lyr - 1, :]
    )
    DS[2 * ice.nbr_lyr - 1, :] = np.nan_to_num(
        E[2 * ice.nbr_lyr - 1, :] / B[2 * ice.nbr_lyr - 1, :]
    )

    # for all layers above bottom layer, starting at second-to-bottom and
    # progressing towards surface:
    # Toon et al Eq 46
    X = np.zeros([ice.nbr_lyr * 2, model_config.nbr_wvl])

    for i in np.arange(2 * ice.nbr_lyr - 2, -1, -1):
        X[i, :] = 1 / (B[i, :] - (D[i, :] * AS[i + 1, :]))
        AS[i, :] = np.nan_to_num(A[i, :] * X[i, :])
        DS[i, :] = np.nan_to_num((E[i, :] - (D[i, :] * DS[i + 1, :])) * X[i, :])

    # then for all layers, progressing from surface to bottom
    # Toon et al Eq 47
    Y = np.zeros([ice.nbr_lyr * 2, model_config.nbr_wvl])

    for i in np.arange(0, 2 * ice.nbr_lyr, 1):
        if i == 0:
            Y[0, :] = DS[0, :]
        else:
            Y[i, :] = DS[i, :] - (AS[i, :] * Y[i - 1, :])

    return Y


def layer_fluxes(
    ice,
    illumination,
    model_config,
    tau_clm,
    tau_star,
    lam,
    Y,
    GAMMA,
    C_pls_btm,
    C_mns_btm,
    C_pls_top,
    C_mns_top,
    e1,
    e2,
    e3,
    e4,
    mu_one,
):
    """Calculate total fluxes.

    Args:
        ice: instance of Ice class
        illumination: instance of Illumination class
        model_config: instance of ModelConfig class
        tau_clm: optical thickness of column
        lam: intermediate coefficient
        GAMMA: intermediate coefficient
        tau_star: delta scaled optical thickness
        C_pls_btm: upwards flux from bottom of layer
        C_mns_btm: downwards flux from bottom of layer
        C_pls_top: upwards flux from top of layer
        C_mns_top: downwards flux from top of layer
        e1: coefficient for matrix solution
        e2: coefficient for matrix solution
        e3: coefficient for matrix solution
        e4: coefficient for matrix solution
        mu_one: adjusted indicence angle

    Returns:
        F_net: net flux in each layer
        F_top_pls: net upwards flux from upper surface of each layer
        F_btm_net: net flux at bottom surface
        F_top_net: net fluxes at upper surface
        intensity: mean intensity at base of each layer

    """
    direct = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    F_net = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    F_btm_net = np.zeros([1, model_config.nbr_wvl])
    F_top_net = np.zeros([1, model_config.nbr_wvl])
    intensity = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    F_top_pls = np.zeros([1, model_config.nbr_wvl])
    F_up = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    F_down = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])

    # ----------------------------------------------------------------------------------
    # CALCULATE DIRECT BEAM FLUX AT BOTTOM OF EACH LAYER

    # loop through layers
    for i in np.arange(0, ice.nbr_lyr, 1):

        # (Toon et al. eq 50)
        direct[i, :] = (
            illumination.mu_not
            * np.pi
            * illumination.Fs
            * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / illumination.mu_not)
        )

        # net flux (positive upward = F_up - F_down) at the base of each layer
        # (Toon et al. Eq 48)
        F_net[i, :] = (
            (Y[2 * i, :] * (e1[i, :] - e3[i, :]))
            + (Y[2 * i + 1, :] * (e2[i, :] - e4[i, :]))
            + C_pls_btm[i, :]
            - C_mns_btm[i, :]
            - direct[i, :]
        )

        # mean intensity at the base of each layer (Toon et al. Eq 49)
        intensity[i, :] = (1 / mu_one) * (
            Y[2 * i, :] * (e1[i, :] + e3[i, :])
            + Y[2 * i + 1, :] * (e2[i, :] + e4[i, :])
            + C_pls_btm[i, :]
            + C_mns_btm[i, :]
        ) + (direct[i, :] / illumination.mu_not)
        intensity[i, :] = intensity[i, :] / (4 * np.pi)

    # Upward flux at upper model boundary (Toon et al Eq 31)
    F_top_pls = (
        (Y[0, :] * (np.exp(-lam[0, :] * tau_star[0, :]) + GAMMA[0, :]))
        + (Y[1, :] * (np.exp(-lam[0, :] * tau_star[0, :]) - GAMMA[0, :]))
        + C_pls_top[0, :]
    )

    for i in np.arange(0, ice.nbr_lyr, 1):

        # Upward flux at the bottom of each layer interface (Toon et al. Eq31)
        F_up[i, :] = (
            Y[2 * i, :]
            * (np.exp(0) + GAMMA[i, :] * np.exp(-lam[i, :] * tau_star[i, :]))
            + Y[2 * i + 1, :]
            * (np.exp(0) - GAMMA[i, :] * np.exp(-lam[i, :] * tau_star[i, :]))
            + C_pls_btm[i, :]
        )

        # Downward flux at the bottom of each layer interface (Toon et al. Eq32)
        # plus direct beam component
        F_down[i, :] = (
            Y[2 * i, :]
            * (GAMMA[i, :] * np.exp(0) + np.exp(-lam[i, :] * tau_star[i, :]))
            + Y[2 * i + 1, :]
            * (GAMMA[i, :] * np.exp(0) - np.exp(-lam[i, :] * tau_star[i, :]))
            + C_mns_btm[i, :]
            + direct[i, :]
        )

    # Net flux at lower model boundary = bulk transmission through entire media
    # = energy absorbed by underlying surface
    F_btm_net[0, :] = -F_net[ice.nbr_lyr - 1, :]

    return F_net, F_top_pls, F_btm_net, F_top_net, intensity


def absorbed_fluxes(ice, model_config, F_net, F_top_net):
    """Calculates energy absorbed in each layer.

    Args:
        ice: instance of Ice class
        model_config: instance of ModelConfig class
        F_net: net flux in each layer
        F_top_net: net flux at upper model boundary

    Returns:
        F_abs: absorbed flux in each layer

    """

    F_abs = np.zeros([ice.nbr_lyr, model_config.nbr_wvl])
    # absorbed flux in each layer (negative if there is net emission (bnd_typ = 4))
    for i in np.arange(0, ice.nbr_lyr, 1):
        if i == 0:
            F_abs[0, :] = F_net[0, :] - F_top_net
        else:
            F_abs[i, :] = F_net[i, :] - F_net[i - 1, :]

    return F_abs


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
    energy_sum = (
        (illumination.mu_not * np.pi * illumination.Fs)
        + illumination.Fd
        - (sum(F_abs) + F_btm_net + F_top_pls)
    )

    # spectrally-integrated terms:
    # energy conservation total error
    energy_error = abs(np.sum(energy_sum))

    if energy_error > 1e-10:
        energy_conservation_error = np.sum(abs(energy_sum))
        raise ValueError(f"CONSERVATION OF ENERGY ERROR OF {energy_conservation_error}")

    return


def get_outputs(
    model_config, ice, illumination, F_top_pls, F_top_net, F_btm_net, F_abs, L_snw
):
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
    intensity_top = np.zeros(model_config.nbr_wvl)
    abs_vis = np.zeros(ice.nbr_lyr)
    abs_nir = np.zeros(ice.nbr_lyr)
    outputs = Outputs()

    # surface planar intensity
    intensity_top[:] = F_top_pls + (
        (illumination.mu_not * np.pi * illumination.Fs) + illumination.Fd
    )

    # Hemispheric wavelength-dependent albedo
    albedo = F_top_pls / (
        (illumination.mu_not * np.pi * illumination.Fs) + illumination.Fd
    )

    # Net flux at upper model boundary
    F_top_net[0, :] = F_top_pls - (
        (illumination.mu_not * np.pi * illumination.Fs) + illumination.Fd
    )

    for i in np.arange(0, ice.nbr_lyr, 1):
        abs_vis[i] = np.sum(F_abs[i, 0 : model_config.vis_max_idx])
        abs_nir[i] = np.sum(
            F_abs[i, model_config.vis_max_idx : model_config.nir_max_idx]
        )

    # Spectrally-integrated absorption in each layer:
    outputs.abs_slr = np.sum(F_abs, axis=1)

    # energy absorbed by underlying substrate
    outputs.abs_slr_btm = sum(np.squeeze(F_btm_net))
    outputs.abs_vis_btm = sum(np.squeeze(F_btm_net[0 : model_config.vis_max_idx]))
    outputs.abs_nir_btm = sum(
        np.squeeze(F_btm_net[0, model_config.vis_max_idx : model_config.nir_max_idx])
    )

    # Calculate radiative heating rate in kelvin per second.
    # Multiply by 3600 to convert to K per hour
    # specfic heta capacity of ice = 2117 J kg-1 K-1
    heat_rt = outputs.abs_slr / (L_snw * 2117)  # [K / s]
    outputs.heat_rt = heat_rt * 3600  # [K / hr]

    # ----------------------------------------------------------------------------------
    # Re-alias results for outputting
    # ----------------------------------------------------------------------------------

    # total incident insolation(Wm - 2)
    outputs.total_insolation = np.sum(
        (illumination.mu_not * np.pi * illumination.Fs) + illumination.Fd
    )
    outputs.albedo = albedo

    # energy absorbed by all snow layers
    outputs.abs_slr_tot = np.sum(np.sum(F_abs))

    # Spectrally - integrated solar, visible, and NIR albedos:
    outputs.BBA = np.sum(illumination.flx_slr * albedo) / np.sum(illumination.flx_slr)

    outputs.BBAVIS = sum(
        illumination.flx_slr[0 : model_config.vis_max_idx]
        * albedo[0 : model_config.vis_max_idx]
    ) / sum(illumination.flx_slr[0 : model_config.vis_max_idx])

    outputs.BBANIR = sum(
        illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx]
        * albedo[model_config.vis_max_idx : model_config.nir_max_idx]
    ) / sum(illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx])

    # % Spectrally - integrated VIS and NIR total snowpack absorption:
    outputs.abs_vis_tot = sum(
        illumination.flx_slr[0 : model_config.vis_max_idx]
        * (1 - albedo[0 : model_config.vis_max_idx])
    )
    outputs.abs_nir_tot = sum(
        illumination.flx_slr[model_config.vis_max_idx : model_config.nir_max_idx]
        * (1 - albedo[model_config.vis_max_idx : model_config.nir_max_idx])
    )
    outputs.absorbed_flux_per_layer = F_abs

    return outputs


def apply_smoothing_function(albedo, model_config):

    yhat = savgol_filter(albedo, model_config.window_size, model_config.poly_order)
    albedo = yhat

    return albedo


if __name__ == "__main__":
    pass
