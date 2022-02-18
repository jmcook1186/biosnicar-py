from setup_snicar import *
from classes import *
import numpy as np


def adding_doubling_solver(
    tau, ssa, g, L_snw, ice, illumination, model_config, rt_config
):

    """
    This script is one of the two optional radiativ transfer solvers available in this
    package. This script deals with the adding-doubling method as translated from
    MATLAB code from Chloe Whicker (UMich) - October 2020. When it becomes available,
    any use of this adding-doubling script should cite Chloe's paper.

    This is the appropriate solver for any configuration where solid ice layers and
    fresnel reflection
    are included.

    """

    outputs = Outputs()

    # DEFINE CONSTANTS AND SET UP ARRAYS

    tau0 = tau.T  # read and transpose tau
    g0 = g.T  # read and transpose g
    ssa0 = ssa.T  # read and transpose ssa

    epsilon = 1e-5  # to deal with singularity
    exp_min = (
        1e-5  # exp(-500)  # minimum number that is not zero - zero will raise error
    )
    trmin = 1e-4  # minimum transmissivity
    puny = 1e-10  # not sure how should we define this

    gauspt = [
        0.9894009,
        0.9445750,
        0.8656312,
        0.7554044,
        0.6178762,
        0.4580168,
        0.2816036,
        0.0950125,
    ]  # gaussian angles (radians)
    gauswt = [
        0.0271525,
        0.0622535,
        0.0951585,
        0.1246290,
        0.1495960,
        0.1691565,
        0.1826034,
        0.1894506,
    ]  # gaussian weights

    vis_max_idx = 50  # index of maximum visible wavelength (0.7 um)
    nir_max_idx = 480  # index of max nir wavelength (5 um)

    # empty arrays
    trndir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trntdr = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trndif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rupdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rupdif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    trndir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rdndif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdirup = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdifup = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdirdn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    fdifdn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    dfdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    dfdif = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rdir = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rdif_a = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    rdif_b = np.zeros(
        shape=[model_config.nbr_wvl, ice.nbr_lyr + 1]
    )  # layer reflectivity to diffuse radiation from below
    tdir = np.zeros(
        shape=[model_config.nbr_wvl, ice.nbr_lyr + 1]
    )  # layer transmission to direct radiation (solar beam + diffuse)
    tdif_a = np.zeros(
        shape=[model_config.nbr_wvl, ice.nbr_lyr + 1]
    )  # layer transmission to diffuse radiation from above
    tdif_b = np.zeros(
        shape=[model_config.nbr_wvl, ice.nbr_lyr + 1]
    )  # layer transmission to diffuse radiation from below
    trnlay = np.zeros(
        shape=[model_config.nbr_wvl, ice.nbr_lyr + 1]
    )  # solar beam transm for layer (direct beam only)
    F_up = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    F_dwn = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr + 1])
    F_abs = np.zeros(shape=[model_config.nbr_wvl, ice.nbr_lyr])
    F_abs_vis = np.zeros(shape=[ice.nbr_lyr])
    F_abs_nir = np.zeros(shape=[ice.nbr_lyr])
    trndir[:, 0] = 1
    trntdr[:, 0] = 1
    trndif[:, 0] = 1
    rdndif[:, 0] = 0
    n_real = np.zeros(shape=480)

    # set the underlying ground albedo
    rupdir[
        :, ice.nbr_lyr
    ] = ice.sfc  # reflectivity to direct radiation for layers below
    rupdif[
        :, ice.nbr_lyr
    ] = ice.sfc  # reflectivity to diffuse radiation for layers below

    # if there are non zeros in layer type, grab the index of the first fresnel layer
    # if there are non-zeros in layer type, load in the precalculated diffuse fresnel
    # reflection
    # (precalculated as large no. of gaussian points required for convergence)
    if np.sum(ice.layer_type) > 0:

        lyrfrsnl = ice.layer_type.index(1)

    else:

        lyrfrsnl = 9999999

    mu0 = illumination.mu_not * np.ones(
        480
    )  # cosine of beam angle is equal to incident beam
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
    n_real = (np.sqrt(2) / 2) * (
        temp1 + (temp2**2 + 4 * ice.ref_idx_re**2 * ice.ref_idx_im**2) ** 0.5
    ) ** 0.5
    nr = n_real
    # . Eq. 20: Briegleb and Light 2007: adjusts beam angle
    # (i.e. this is Snell's Law for refraction at interface between media)
    # mu0n = -1 represents light travelling vertically upwards and mu0n = +1
    # represents light travellign vertically downwards
    # mu0n = np.sqrt(1-((1-mu0**2)/(ref_indx*ref_indx)))  (original, before update
    # for diffuse Fresnel reflection)
    mu0n = np.cos(
        np.arcsin(np.sin(np.arccos(mu0)) / nr)
    )  # this version accounts for diffuse fresnel reflection

    # proceed down one layer at a time: if the total transmission to
    # the interface just above a given layer is less than trmin, then no
    # Delta-Eddington computation for that layer is done.

    for lyr in np.arange(0, ice.nbr_lyr, 1):  # loop through layers

        # condition: if current layer is above fresnel layer or the
        # top layer is a Fresnel layer
        if lyr < lyrfrsnl or lyrfrsnl == 0:

            mu0n = mu0

        else:
            # within or below fl
            mu0n = mu0n

        # calculation over layers with penetrating radiation
        # includes optical thickness, single scattering albedo,
        # asymmetry parameter and total flux
        tautot = tau0[:, lyr]
        wtot = ssa0[:, lyr]
        gtot = g0[:, lyr]
        ftot = g0[:, lyr] * g0[:, lyr]

        # coefficient for delta eddington solution for all layers
        # Eq. 50: Briegleb and Light 2007
        ts = (1 - (wtot * ftot)) * tautot  # layer delta-scaled extinction optical depth
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
        tdif_a[:, lyr] = 4 * ue / ne  # T BAR layer transmissivity to DIFFUSE radiation

        # evaluate rdir, tdir for direct beam
        trnlay[:, lyr] = np.maximum(
            np.full((model_config.nbr_wvl,), exp_min), np.exp(-ts / mu0n)
        )  # transmission from TOA to interface

        #  Eq. 50: Briegleb and Light 2007  alpha and gamma for direct radiation
        alp = (
            (0.75 * ws * mu0n)
            * (1 + gs * (1 - ws))
            / (1 - lm**2 * mu0n**2 + epsilon)
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

        # recalculate rdif,tdif using direct angular integration over rdir,tdir,
        # since Delta-Eddington rdif formula is not well-behaved (it is usually
        # biased low and can even be negative)  use ngmax angles and gaussian
        # integration for most accuracy:

        R1 = rdif_a[:, lyr]  # use R1 as temporary var
        T1 = tdif_a[:, lyr]  # use T1 as temporary var
        swt = 0
        smr = 0
        smt = 0

        # loop through the gaussian angles for the AD integral
        for ng in np.arange(0, len(gauspt), 1):  # gaussian angles (radians)

            mu = gauspt[ng]  # solar zenith angles
            gwt = gauswt[ng]  # gaussian weight
            swt = swt + mu * gwt  # sum of weights
            trn = np.maximum(
                np.full((model_config.nbr_wvl,), exp_min), np.exp(-ts / mu)
            )  # transmission

            alp = (
                (0.75 * ws * mu)
                * (1 + gs * (1 - ws))
                / (1 - lm**2 * mu**2 + epsilon)
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

        rdif_a[:, lyr] = smr / swt
        tdif_a[:, lyr] = smt / swt

        # homogeneous layer
        rdif_b[:, lyr] = rdif_a[:, lyr]
        tdif_b[:, lyr] = tdif_a[:, lyr]

        # ------------------------------------------------------------------------------
        # Fresnel layer
        # ------------------------------------------------------------------------------

        if lyr == lyrfrsnl:
            ref_indx = (
                ice.ref_idx_re + 1j * ice.ref_idx_im
            )  # np.complex(ice.ref_idx_re,ice.ref_idx_im)
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

                Rf_dif_a = ice.fl_r_dif_a[
                    wl
                ]  # reflection from diffuse unpolarized radiation
                Tf_dif_a = (
                    1 - Rf_dif_a
                )  # transmission from diffuse unpolarized radiation

                # diffuse reflection of flux arriving from below
                Rf_dif_b = ice.fl_r_dif_b[wl]
                Tf_dif_b = 1 - Rf_dif_b

                # -----------------------------------------------------------------------
                # the lyr = lyrfrsnl layer properties are updated to combine
                # the fresnel (refractive) layer, always taken to be above
                # the present layer lyr (i.e. be the top interface):

                rintfc = 1 / (
                    1 - Rf_dif_b * rdif_a[wl, lyr]
                )  # denom interface scattering

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
                rdif_a[wl, lyr] = (
                    Rf_dif_a + Tf_dif_a * rdif_a[wl, lyr] * rintfc * Tf_dif_b
                )

                # R BAR = layer reflectivity to DIFFUSE radiation (below)
                # Eq. B10  Briegleb & Light 2007
                rdif_b[wl, lyr] = (
                    rdif_b[wl, lyr]
                    + tdif_b[wl, lyr] * Rf_dif_b * rintfc * tdif_a[wl, lyr]
                )

                # T BAR layer transmissivity to DIFFUSE radiation (above),
                # Eq. B9  Briegleb & Light 2007
                tdif_a[wl, lyr] = tdif_a[wl, lyr] * rintfc * Tf_dif_a

                # Eq. B10  Briegleb & Light 2007
                tdif_b[wl, lyr] = tdif_b[wl, lyr] * rintfc * Tf_dif_b

                # update trnlay to include fresnel transmission
                trnlay[wl, lyr] = Tf_dir_a * trnlay[wl, lyr]

                # end lyr = lyrfrsnl condition
                # end trntdr[lyr, wl] > trmin condition

            #   Calculate the solar beam transmission, total transmission, and
            #   reflectivity for diffuse radiation from below at interface lyr,
            #   the top of the current layer lyr:
            #
            #                layers       interface
            #
            #         ---------------------  lyr-1
            #                  lyr-1
            #         ---------------------  lyr
            #                   lyr
            #         ---------------------
            #   note that we ignore refraction between sea ice and underlying ocean:
            #
            #                layers       interface
            #
            #         ---------------------  lyr-1
            #                  lyr-1
            #         ---------------------  lyr
            #         \\\\\\\ ocean \\\\\\\

        # Eq. 51  Briegleb and Light 2007

        trndir[:, lyr + 1] = (
            trndir[:, lyr] * trnlay[:, lyr]
        )  # solar beam transmission from top
        # trnlay = exp(-ts/illumination.mu_not) = direct solar beam transmission

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
        rdndif[:, lyr + 1] = rdif_b[:, lyr] + (
            tdif_b[:, lyr] * rdndif[:, lyr] * refkm1 * tdif_a[:, lyr]
        )  # reflectivity to diffuse radiation for layers above
        trndif[:, lyr + 1] = (
            trndif[:, lyr] * refkm1 * tdif_a[:, lyr]
        )  # diffuse transmission to diffuse beam for layers above

    # end main level loop  number of layers

    # ! compute reflectivity to direct and diffuse radiation for layers
    # ! below by adding succesive layers starting from the underlying
    # ! ocean and working upwards:
    # !
    # !              layers       interface
    # !
    # !       ---------------------  lyr
    # !                 lyr
    # !       ---------------------  lyr+1
    # !                lyr+1
    # !       ---------------------

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

    # fluxes at interface

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

        if np.min(dfdir[:, lyr]) < puny:

            dfdir[:, lyr] = np.zeros(
                (model_config.nbr_wvl,), dtype=int
            )  # echmod necessary?
            # dfdif = fdifdn - fdifup

        dfdif[:, lyr] = trndif[:, lyr] * (1 - rupdif[:, lyr]) * refk

        if np.min(dfdif[:, lyr]) < puny:

            dfdif[:, lyr] = np.zeros(
                (model_config.nbr_wvl,), dtype=int
            )  #!echmod necessary?

    # ----- End Radiative Solver Adding Doubling Method -----

    # ----- Calculate fluxes ----

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

    # Absorbed flux in each layer
    F_abs[:, :] = F_net[:, 1:] - F_net[:, :-1]

    # # albedo
    # acal  = F_up[:,0]/F_dwn[:,0]

    # Upward flux at upper model boundary
    F_top_pls = F_up[:, 0]

    # Net flux at lower model boundary = bulk transmission through entire
    # media = absorbed radiation by underlying surface:
    F_btm_net = -F_net[:, ice.nbr_lyr]

    # Spectrally-integrated absorption in each layer:
    F_abs_slr = np.sum(F_abs, axis=0)

    for i in np.arange(0, ice.nbr_lyr, 1):  # [0,1,2,3,4]

        F_abs_vis[i] = sum(F_abs[0:vis_max_idx, i])
        F_abs_nir[i] = sum(F_abs[vis_max_idx:nir_max_idx, i])

    # Spectrally-integrated absorption by underlying surface:
    outputs.F_abs_btm = np.sum(F_btm_net, axis=0)
    outputs.F_abs_vis_btm = np.sum(F_btm_net[0:vis_max_idx], axis=0)
    outputs.F_abs_nir_btm = np.sum(F_btm_net[vis_max_idx : nir_max_idx + 1], axis=0)

    # Radiative heating rate:
    heat_rt = F_abs_slr / (L_snw * 2117)  # [K/s] 2117 = specific heat ice (J kg-1 K-1)
    outputs.heat_rt = heat_rt * 3600  # [K/hr]

    # Energy conservation check:
    # Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
    energy_sum = (
        (illumination.mu_not * np.pi * illumination.Fs)
        + illumination.Fd
        - (np.sum(F_abs, axis=1) + F_btm_net + F_top_pls)
    )

    energy_conservation_error = sum(abs(energy_sum))

    if energy_conservation_error > 1e-10:

        raise ValueError(f"energy conservation error: {energy_conservation_error}")

    albedo = F_up[:, 0] / F_dwn[:, 0]

    # Spectrally-integrated solar, visible, and NIR albedos:

    alb_bb = np.sum(illumination.flx_slr * albedo) / np.sum(illumination.flx_slr)

    outputs.BBAVIS = np.sum(
        illumination.flx_slr[0:vis_max_idx] * albedo[0:vis_max_idx]
    ) / np.sum(illumination.flx_slr[0:vis_max_idx])

    outputs.BBANIR = np.sum(
        illumination.flx_slr[vis_max_idx:nir_max_idx] * albedo[vis_max_idx:nir_max_idx]
    ) / np.sum(illumination.flx_slr[vis_max_idx:nir_max_idx])

    # Spectrally-integrated VIS and NIR total snowpack absorption:
    outputs.abs_vis_tot = np.sum(illumination.flx_slr[0:vis_max_idx] * (1 - albedo[0:vis_max_idx]))
    outputs.abs_nir_tot = np.sum(
        illumination.flx_slr[vis_max_idx:nir_max_idx]
        * (1 - albedo[vis_max_idx:nir_max_idx])
    )

    # -------------------- OUTPUT  --------------------

    flx_dwn_spc = (
        illumination.mu_not * np.pi * illumination.Fs + illumination.Fd
    )  # spectral downwelling flux at model top [W/m2/band]
    outputs.BBA = alb_bb  # solar broadband albedo
    outputs.abs_slr_tot = np.sum(F_abs_slr)  # total solar absorption by entire snow column
    #       (not including underlying substrate) [W/m2]
    # abs_snw_vis = np.sum(F_abs_vis)
    #   visible solar absorption by entire snow column
    #   (not including underlying substrate) [W/m2]
    # abs_snw_nir = np.sum(F_abs_nir)
    #   near-IR solar absorption by entire snow column
    #   (not including underlying substrate) [W/m2]
    # abs_spc = np.sum(F_abs,axis=1)
    #  spectral absorption by entire snow column [W/m2/band]

    # abs_snw_top_slr = F_abs_slr[0]        # top snow layer solar absorption [W/m2]
    # abs_snw_top_vis = F_abs_vis[0]        # top snow layer VIS absorption [W/m2]
    # abs_snw_top_nir = F_abs_nir[0]        # top snow layer NIR absorption [W/m2]

    # abs_ground_slr  = F_abs_btm
    #    total solar absorption by underlying substrate [W/m2]
    # abs_ground_vis  = F_abs_vis_btm
    #    visible absorption by underlying substrate [W/m2]
    # abs_ground_nir  = F_abs_nir_btm
    #   near-IR absorption by underlying substrate [W/m2]
    outputs.albedo = albedo

    return outputs
