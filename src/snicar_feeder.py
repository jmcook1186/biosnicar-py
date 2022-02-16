



    # ----------------------------------------------------------------------------------
    # Read in impurity optical properties
    # ----------------------------------------------------------------------------------

    # Load optical properties ssa, MAC and g
    # (one row per impurity, one column per wvalengths)
    # Load mass concentrations MSS per layer
    # (one row per layer, one column per impurity)

    ssa_aer = np.zeros([Inputs.nbr_aer, nbr_wvl])
    mac_aer = np.zeros([Inputs.nbr_aer, nbr_wvl])
    g_aer = np.zeros([Inputs.nbr_aer, nbr_wvl])
    mss_aer = np.zeros([Inputs.nbr_lyr, Inputs.nbr_aer])

    for aer in range(Inputs.nbr_aer):

        impurity_properties = xr.open_dataset(str(dir_lap_files + files[aer]))

        g_aer[aer, :] = impurity_properties["asm_prm"].values
        ssa_aer[aer, :] = impurity_properties["ss_alb"].values

        # coated particles: use ext_cff_mss_ncl for MAC
        if files[aer] == Inputs.file_brwnC2 or files[aer] == Inputs.file_soot2:
            mac_aer[aer, :] = impurity_properties["ext_cff_mss_ncl"].values

        else:
            mac_aer[aer, :] = impurity_properties["ext_cff_mss"].values

        if files[aer] == Inputs.file_glacier_algae:
            # if GA_units == 1, GA concentration provided in cells/mL
            # mss_aer should be in cells/kg
            # thus mss_aer is divided by kg/mL ice = 917*10**(-6)
            # with density of ice 917 kg m3
            if Inputs.GA_units == 1:

                mss_aer[0:Inputs.nbr_lyr, aer] = np.array(mass_concentrations[aer]) / (
                    917 * 10 ** (-6)
                )

            else:
                mss_aer[0:Inputs.nbr_lyr, aer] = np.array(mass_concentrations[aer]) * 1e-9

        elif files[aer] == Inputs.file_snw_alg:
            # if SA_units == 1, SA concentration provided in cells/mL
            # but mss_aer should be in cells/kg
            # thus mss_aer is divided by kg/mL ice = 917*10**(-6)
            # with density of ice 917 kg m3
            if Inputs.SA_units == 1:
                mss_aer[0:Inputs.nbr_lyr, aer] = np.array(mass_concentrations[aer]) / (
                    917 * 10 ** (-6)
                )

            else:
                mss_aer[0:Inputs.nbr_lyr, aer] = np.array(mass_concentrations[aer]) * 1e-9

        else:
            # conversion to kg/kg ice from ng/g
            mss_aer[0:Inputs.nbr_lyr, aer] = np.array(mass_concentrations[aer]) * 1e-9

        # if c_factor provided, then mss_aer multiplied by c_factor
        if (
            files[aer] == Inputs.file_glacier_algae
            and isinstance(Inputs.c_factor_GA, (int, float))
            and (Inputs.c_factor_GA > 0)
        ):
            mss_aer[0:Inputs.nbr_lyr, aer] = Inputs.c_factor_GA * mss_aer[0:Inputs.nbr_lyr, aer]

        if (
            files[aer] == Inputs.file_snw_alg
            and isinstance(Inputs.c_factor_SA, (int, float))
            and (Inputs.c_factor_SA > 0)
        ):
            mss_aer[0:Inputs.nbr_lyr, aer] = Inputs.c_factor_SA * mss_aer[0:Inputs.nbr_lyr, aer]

    # ----------------------------------------------------------------------------------
    # Begin solving Radiative Transfer
    # -----------------------------------------------------------------------------------

    # 1. Calculate effective tau (optical Inputs.depth),
    # ssa (single scattering albedo) and
    # g (assymetry parameter) for the ice +
    # impurities mixture.

    # ssa and g for the individual components has
    # been calculated using Mie theory and
    # stored in a netcdf file. Here, these values
    # are combined to give an overall
    # ssa and g for the ice + impurity mixture

    # initialize arrays
    g_sum = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    ssa_sum = np.zeros([Inputs.nbr_lyr, Inputs.nbr_aer, nbr_wvl])
    tau = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    ssa = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    g = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    L_aer = np.zeros([Inputs.nbr_lyr, Inputs.nbr_aer])
    tau_aer = np.zeros([Inputs.nbr_lyr, Inputs.nbr_aer, nbr_wvl])
    tau_sum = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    ssa_sum = np.zeros([Inputs.nbr_lyr, nbr_wvl])
    L_snw = np.zeros(Inputs.nbr_lyr)
    tau_snw = np.zeros([Inputs.nbr_lyr, nbr_wvl])

    # for each layer, the layer mass (L) is density * layer thickness
    # for each layer the optical Inputs.depth is
    # the layer mass * the mass extinction coefficient
    # first for the ice in each layer

    for i in range(Inputs.nbr_lyr):

        L_snw[i] = Inputs.rho_layers[i] * Inputs.dz[i]

        for j in range(Inputs.nbr_aer):

            # kg ice m-2 * cells kg-1 ice = cells m-2
            L_aer[i, j] = L_snw[i] * mss_aer[i, j]
            # cells m-2 * m2 cells-1
            tau_aer[i, j, :] = L_aer[i, j] * mac_aer[j, :]
            tau_sum[i, :] = tau_sum[i, :] + tau_aer[i, j, :]
            ssa_sum[i, :] = ssa_sum[i, :] + (tau_aer[i, j, :] * ssa_aer[j, :])
            g_sum[i, :] = g_sum[i, :] + (tau_aer[i, j, :] * ssa_aer[j, :] * g_aer[j, :])

            # ice mass = snow mass - impurity mass (generally tiny correction)
            # if aer == algae and L_aer is in cells m-2, should be converted
            # to m-2 kg-1 : 1 cell = 1ng = 10**(-12) kg

            if files[j] == Inputs.file_glacier_algae and Inputs.GA_units == 1:

                L_snw[i] = L_snw[i] - L_aer[i, j] * 10 ** (-12)

            elif files[j] == Inputs.file_snw_alg and Inputs.SA_units == 1:

                L_snw[i] = L_snw[i] - L_aer[i, j] * 10 ** (-12)

            else:

                L_snw[i] = L_snw[i] - L_aer[i, j]

        tau_snw[i, :] = L_snw[i] * MAC_snw[i, :]
        # finally, for each layer calculate the effective ssa, tau and g
        # for the snow+LAP
        tau[i, :] = tau_sum[i, :] + tau_snw[i, :]
        ssa[i, :] = (1 / tau[i, :]) * (ssa_sum[i, :] + (ssa_snw[i, :] * tau_snw[i, :]))
        g[i, :] = (1 / (tau[i, :] * (ssa[i, :]))) * (
            g_sum[i, :] + (g_snw[i, :] * ssa_snw[i, :] * tau_snw[i, :])
        )

    Inputs.tau = tau
    Inputs.ssa = ssa
    Inputs.g = g
    Inputs.L_snw = L_snw

    # just in case any unrealistic values arise (none detected so far)
    ssa[ssa <= 0] = 0.00000001
    ssa[ssa >= 1] = 0.99999999
    g[g <= 0] = 0.00001
    g[g >= 1] = 0.99999

    Inputs.fl_r_dif_a = fl_r_dif_a
    Inputs.fl_r_dif_b = fl_r_dif_b
    Inputs.refidx_re = refidx_re
    Inputs.refidx_im = refidx_im
    Inputs.flx_slr = flx_slr
    Inputs.Fs = Fs
    Inputs.Fd = Fd
    Inputs.tau = tau
    Inputs.ssa = ssa
    Inputs.g = g
    Inputs.L_snw = L_snw

    # CALL RT SOLVER (toon_solver  = toon ET AL, TRIDIAGONAL MATRIX METHOD;
    # add_double = ADDING-DOUBLING METHOD)

    Outputs = c.namedtuple(
        "Outputs",
        ["wvl", "albedo", "BBA", "BBAVIS", "BBANIR", "abs_slr", "heat_rt", "abs_ice"],
    )

    if Inputs.toon:

        (
            Outputs.wvl,
            Outputs.albedo,
            Outputs.BBA,
            Outputs.BBAVIS,
            Outputs.BBANIR,
            Outputs.abs_slr,
            Outputs.heat_rt,
        ) = toon_solver.toon_solver(Inputs)

    if Inputs.add_double:

        (
            Outputs.wvl,
            Outputs.albedo,
            Outputs.BBA,
            Outputs.BBAVIS,
            Outputs.BBANIR,
            Outputs.abs_slr,
            Outputs.heat_rt,
        ) = adding_doubling.adding_doubling_solver(Inputs)

    return Outputs
