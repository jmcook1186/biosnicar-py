def toon_solver(Inputs):
    """
    Three 2-stream approximations are available: Eddington,
    Quadrature and hemispheric mean. The equations for each
    approximation are provided in Toon et al. (1989) Table 1.

    The hemispheric mean scheme is derived by assuming that the
    phase function is equal to 1  + g  in the forward scattering
    hemisphere and to 1  - g  in the backward scattering hemisphere.
    The asymmetry parameter is g. The hemispheric mean is only
    useful for infrared wavelengths
    """

    # load variables from input table

    import numpy as np

    direct = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    F_net = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    F_btm_net = np.zeros([1, Inputs.nbr_wvl])
    F_top_net = np.zeros([1, Inputs.nbr_wvl])
    intensity = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    F_top_pls = np.zeros([1, Inputs.nbr_wvl])
    F_up = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    F_down = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    F_net2 = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    intensity2 = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    intensity2_top = np.zeros(Inputs.nbr_wvl)
    F_abs = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    abs_vis = np.zeros(Inputs.nbr_lyr)
    abs_nir = np.zeros(Inputs.nbr_lyr)

    # ----------------------------------------------------------------------------------
    # PERFORM DELTA TRANSFORMATION IF REQUIRED
    # ----------------------------------------------------------------------------------
    # The star represents the delta transformed quantity
    # if no delta transformation is applied, the starred quantity
    # is equal to the unstarred quantity

    if Inputs.delta:
        g_star = Inputs.g / (1 + Inputs.g)
        ssa_star = ((1 - (Inputs.g**2)) * Inputs.ssa) / (
            1 - (Inputs.ssa * (Inputs.g**2))
        )
        tau_star = (1 - (Inputs.ssa * (Inputs.g**2))) * Inputs.tau

    else:
        g_star = Inputs.g
        ssa_star = Inputs.ssa
        tau_star = Inputs.tau

    # CALCULATE TOTAL OPTICAL DEPTH OF ENTIRE COLUMN
    # i.e. tau_clm = total optical depth from upper boundary
    # to upper boundary of layer n. This is therefore a cumulative
    # quantity - subsequently lower layers contain the sum of the
    # # optical depth of all overlying layers

    tau_clm = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    for i in np.arange(1, Inputs.nbr_lyr, 1):
        # start loop from 2nd layer, i.e. index = 1
        tau_clm[i, :] = tau_clm[i - 1, :] + tau_star[i - 1, :]

    # SET BOUNDARY CONDITION: BOTTOM BOUNDARY
    # calculate radiation reflected skywards by underlying surface
    # (i.e. lower model boundary)
    # remainder is lost

    S_sfc = (
        Inputs.R_sfc
        * Inputs.mu_not
        * np.exp(
            -(tau_clm[Inputs.nbr_lyr - 1, :] + tau_star[Inputs.nbr_lyr - 1, :])
            / Inputs.mu_not
        )
        * np.pi
        * Inputs.Fs
    )

    # ----------------------------------------------------------------------------------
    # Apply Two-Stream Approximation (Toon et al, table 1)
    # ----------------------------------------------------------------------------------

    if Inputs.aprx_typ == 1:
        # apply Eddington approximation
        gamma1 = (7 - (ssa_star * (4 + (3 * g_star)))) / 4
        gamma2 = -(1 - (ssa_star * (4 - (3 * g_star)))) / 4
        gamma3 = (2 - (3 * g_star * Inputs.mu_not)) / 4
        gamma4 = 1 - gamma3
        mu_one = 0.5

    elif Inputs.aprx_typ == 2:
        # apply quadrature approximation
        gamma1 = np.sqrt(3) * (2 - (ssa_star * (1 + g_star))) / 2
        gamma2 = ssa_star * np.sqrt(3) * (1 - g_star) / 2
        gamma3 = (1 - (np.sqrt(3) * g_star * Inputs.mu_not)) / 2
        gamma4 = 1 - gamma3
        mu_one = 1 / np.sqrt(3)

    elif Inputs.aprx_typ == 3:
        # apply hemispheric mean approximation
        gamma1 = 2 - (ssa_star * (1 + g_star))
        gamma2 = ssa_star * (1 - g_star)
        gamma3 = (1 - (np.sqrt(3) * g_star * Inputs.mu_not)) / 2
        gamma4 = 1 - gamma3
        mu_one = 0.5

    else:
        raise ValueError("INVALID APRX_TYP: CHOOSE 1, 2 or 3")

    # Toon et al equation 21 and 22
    # Note that the values of lam and GAMMA depend upon gamma1 and gamma2, which
    # vary depending upon the two-stream approximation used
    # variable "lambda" renamed "lam" to avoid confusion with lambda function
    lam = np.sqrt(abs((gamma1**2) - (gamma2**2)))
    GAMMA = gamma2 / (gamma1 + lam)

    # calculate coefficients required for tridiagonal matrix calculation
    # (Toon et al Equation 44)
    e1 = 1 + (GAMMA * np.exp(-lam * tau_star))
    e2 = 1 - (GAMMA * np.exp(-lam * tau_star))
    e3 = GAMMA + np.exp(-lam * tau_star)
    e4 = GAMMA - np.exp(-lam * tau_star)

    # ----------------------------------------------------------------------------------
    # Calculate C-functions
    # ----------------------------------------------------------------------------------

    # C is the direct beam flux calculated at the top and bottom of each layer, i,
    # see Toon equations 23 and 24

    """ N.B. consider adding in stability check here as per Flanner's Matlab code """

    C_pls_btm = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    C_mns_btm = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    C_pls_top = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])
    C_mns_top = np.zeros([Inputs.nbr_lyr, Inputs.nbr_wvl])

    for i in np.arange(0, Inputs.nbr_lyr, 1):

        if np.sum(Inputs.Fs) > 0.0:

            C_pls_btm[i, :] = (
                ssa_star[i, :]
                * np.pi
                * Inputs.Fs
                * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / Inputs.mu_not)
                * (
                    ((gamma1[i, :] - (1 / Inputs.mu_not)) * gamma3[i, :])
                    + (gamma4[i, :] * gamma2[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / (Inputs.mu_not**2)))

            C_mns_btm[i, :] = (
                ssa_star[i, :]
                * np.pi
                * Inputs.Fs
                * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / Inputs.mu_not)
                * (
                    ((gamma1[i, :] + (1 / Inputs.mu_not)) * gamma4[i, :])
                    + (gamma2[i, :] * gamma3[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / Inputs.mu_not**2))

            C_pls_top[i, :] = (
                ssa_star[i, :]
                * np.pi
                * Inputs.Fs
                * np.exp(-tau_clm[i, :] / Inputs.mu_not)
                * (
                    (gamma1[i, :] - (1 / Inputs.mu_not)) * gamma3[i, :]
                    + (gamma4[i, :] * gamma2[i, :])
                )
            ) / ((lam[1, :] ** 2) - (1 / Inputs.mu_not**2))

            C_mns_top[i, :] = (
                ssa_star[i, :]
                * np.pi
                * Inputs.Fs
                * np.exp(-tau_clm[i, :] / Inputs.mu_not)
                * (
                    (gamma1[i, :] + (1 / Inputs.mu_not)) * gamma4[i, :]
                    + (gamma2[i, :] * gamma3[i, :])
                )
            ) / ((lam[i, :] ** 2) - (1 / Inputs.mu_not**2))

        else:
            # no direct-beam flux:
            C_pls_btm[i, :] = 0
            C_mns_btm[i, :] = 0
            C_pls_top[i, :] = 0
            C_mns_top[i, :] = 0

    # Toon equations 41-43.
    # Boundary values for i=1 and i=2Inputs.nbr_lyr, specifics for i=odd and i=even
    # Set up lists
    A = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])
    B = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])
    D = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])
    E = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])

    # ----------------------------------------------------------------------------------
    # Initialize tridiagonal matrix solution
    # ----------------------------------------------------------------------------------

    # expanding the number of layers to 2*Inputs.nbr_lyr so that fluxes at upper and lower
    # layer boundaries can be resolved. This section was confusing to code -
    # for each layer
    # index (n) a second pair of indices (2 x i) are required. Different solutions are
    # applied depending upon whether i is even or odd. To translate the indexing for
    # this from FORTRAN/MATLAB into Python, it was necessary to assert n = (i/2)-1
    # for even layers and n = floor(i/2) for odd layers, with specific rules for the
    # boundaries i = 0 and i = Inputs.nbr_lyrs-1 (i.e. top surface and bottom surface).

    for i in np.arange(0, 2 * Inputs.nbr_lyr, 1):

        # TOP LAYER
        if i == 0:
            A[0, :] = 0.0
            B[0, :] = e1[0, :]
            D[0, :] = -e2[0, :]
            E[0, :] = Inputs.Fd - C_mns_top[0, :]

        # BOTTOM LAYER
        elif i == 2 * Inputs.nbr_lyr - 1:
            A[i, :] = e1[Inputs.nbr_lyr - 1, :] - (
                Inputs.R_sfc * e3[Inputs.nbr_lyr - 1, :]
            )
            B[i, :] = e2[Inputs.nbr_lyr - 1, :] - (
                Inputs.R_sfc * e4[Inputs.nbr_lyr - 1, :]
            )
            D[i, :] = 0.0
            E[i, :] = (
                S_sfc[:]
                - C_pls_btm[Inputs.nbr_lyr - 1, :]
                + (Inputs.R_sfc * C_mns_btm[Inputs.nbr_lyr - 1, :])
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
        elif (i % 2 == 1) and (i < 2 * Inputs.nbr_lyr - 1):

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
    AS = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])
    DS = np.zeros([2 * Inputs.nbr_lyr, Inputs.nbr_wvl])

    np.seterr(divide="ignore", invalid="ignore")
    AS[2 * Inputs.nbr_lyr - 1, :] = np.nan_to_num(
        A[2 * Inputs.nbr_lyr - 1, :] / B[2 * Inputs.nbr_lyr - 1, :]
    )
    DS[2 * Inputs.nbr_lyr - 1, :] = np.nan_to_num(
        E[2 * Inputs.nbr_lyr - 1, :] / B[2 * Inputs.nbr_lyr - 1, :]
    )

    # for all layers above bottom layer, starting at second-to-bottom and
    # progressing towards surface:
    # Toon et al Eq 46
    X = np.zeros([Inputs.nbr_lyr * 2, Inputs.nbr_wvl])
    for i in np.arange(2 * Inputs.nbr_lyr - 2, -1, -1):
        X[i, :] = 1 / (B[i, :] - (D[i, :] * AS[i + 1, :]))
        AS[i, :] = np.nan_to_num(A[i, :] * X[i, :])
        DS[i, :] = np.nan_to_num((E[i, :] - (D[i, :] * DS[i + 1, :])) * X[i, :])

    # then for all layers, progressing from surface to bottom
    # Toon et al Eq 47
    Y = np.zeros([Inputs.nbr_lyr * 2, Inputs.nbr_wvl])

    for i in np.arange(0, 2 * Inputs.nbr_lyr, 1):
        if i == 0:
            Y[0, :] = DS[0, :]
        else:
            Y[i, :] = DS[i, :] - (AS[i, :] * Y[i - 1, :])

    # ----------------------------------------------------------------------------------
    # CALCULATE DIRECT BEAM FLUX AT BOTTOM OF EACH LAYER

    # loop through layers
    for i in np.arange(0, Inputs.nbr_lyr, 1):

        # (Toon et al. eq 50)
        direct[i, :] = (
            Inputs.mu_not
            * np.pi
            * Inputs.Fs
            * np.exp(-(tau_clm[i, :] + tau_star[i, :]) / Inputs.mu_not)
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
        ) + (direct[i, :] / Inputs.mu_not)
        intensity[i, :] = intensity[i, :] / (4 * np.pi)

    # Upward flux at upper model boundary (Toon et al Eq 31)
    F_top_pls = (
        (Y[0, :] * (np.exp(-lam[0, :] * tau_star[0, :]) + GAMMA[0, :]))
        + (Y[1, :] * (np.exp(-lam[0, :] * tau_star[0, :]) - GAMMA[0, :]))
        + C_pls_top[0, :]
    )

    for i in np.arange(0, Inputs.nbr_lyr, 1):
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

        # Derived net (upward-downward) flux (should equal F_net)
        F_net2[i, :] = F_up[i, :] - F_down[i, :]

        intensity2[i, :] = F_up[i, :] + F_down[i, :]

    # surface planar intensity
    intensity2_top[:] = F_top_pls + ((Inputs.mu_not * np.pi * Inputs.Fs) + Inputs.Fd)

    # Net flux at lower model boundary = bulk transmission through entire media
    # = energy absorbed by underlying surface
    F_btm_net[0, :] = -F_net[Inputs.nbr_lyr - 1, :]

    # Hemispheric wavelength-dependent albedo
    albedo = F_top_pls / ((Inputs.mu_not * np.pi * Inputs.Fs) + Inputs.Fd)

    # Net flux at upper model boundary
    F_top_net[0, :] = F_top_pls - ((Inputs.mu_not * np.pi * Inputs.Fs) + Inputs.Fd)

    # absorbed flux in each layer (negative if there is net emission (bnd_typ = 4))
    for i in np.arange(0, Inputs.nbr_lyr, 1):
        if i == 0:
            F_abs[0, :] = F_net[0, :] - F_top_net
        else:
            F_abs[i, :] = F_net[i, :] - F_net[i - 1, :]

    # set indices for constraining calculations to VIS and NIR bands
    vis_max_idx = 39
    nir_max_idx = len(Inputs.wvl)

    # Spectrally-integrated absorption in each layer:
    abs_slr = np.sum(F_abs, axis=1)

    for i in np.arange(0, Inputs.nbr_lyr, 1):
        abs_vis[i] = np.sum(F_abs[i, 0:vis_max_idx])
        abs_nir[i] = np.sum(F_abs[i, vis_max_idx:nir_max_idx])

    # Spectrally - integrated absorption by underlying surface:
    abs_slr_btm = sum(np.squeeze(F_btm_net))
    abs_vis_btm = sum(np.squeeze(F_btm_net[0:vis_max_idx]))
    abs_nir_btm = sum(np.squeeze(F_btm_net[0, vis_max_idx:nir_max_idx]))

    # Calculate radiative heating rate in kelvin per second.
    # Multiply by 3600 to convert to K per hour
    # specfic heta capacity of ice = 2117 J kg-1 K-1
    heat_rt = abs_slr / (Inputs.L_snw * 2117)  # [K / s]
    heat_rt = heat_rt * 3600  # [K / hr]

    # Energy conservation check:
    # % Incident direct + diffuse radiation equals(absorbed + transmitted +
    # bulk_reflected)
    energy_sum = (
        (Inputs.mu_not * np.pi * Inputs.Fs)
        + Inputs.Fd
        - (sum(F_abs) + F_btm_net + F_top_pls)
    )

    # spectrally-integrated terms:
    # energy conservation total error
    energy_error = abs(np.sum(energy_sum))

    if energy_error > 1e-10:
        energy_conservation_error = np.sum(abs(energy_sum))
        raise ValueError(f"CONSERVATION OF ENERGY ERROR OF {energy_conservation_error}")

    # ----------------------------------------------------------------------------------
    # Re-alias results for outputting
    # ----------------------------------------------------------------------------------

    # total incident insolation(Wm - 2)
    total_insolation = np.sum((Inputs.mu_not * np.pi * Inputs.Fs) + Inputs.Fd)

    # energy absorbed by all snow layers
    abs_slr_tot = np.sum(np.sum(F_abs))

    # energy absorbed by underlying substrate
    energy_abs_under_sfc = np.sum(F_btm_net)

    # Spectrally - integrated solar, visible, and NIR albedos:
    BBA = np.sum(Inputs.flx_slr * albedo) / np.sum(Inputs.flx_slr)

    BBAVIS = sum(Inputs.flx_slr[0:vis_max_idx] * albedo[0:vis_max_idx]) / sum(
        Inputs.flx_slr[0:vis_max_idx]
    )

    BBANIR = sum(
        Inputs.flx_slr[vis_max_idx:nir_max_idx] * albedo[vis_max_idx:nir_max_idx]
    ) / sum(Inputs.flx_slr[vis_max_idx:nir_max_idx])

    # % Spectrally - integrated VIS and NIR total snowpack absorption:
    abs_vis_tot = sum(Inputs.flx_slr[0:vis_max_idx] * (1 - albedo[0:vis_max_idx]))
    abs_nir_tot = sum(
        Inputs.flx_slr[vis_max_idx:nir_max_idx] * (1 - albedo[vis_max_idx:nir_max_idx])
    )

    return Inputs.wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, heat_rt
