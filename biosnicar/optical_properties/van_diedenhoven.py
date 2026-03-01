"""Van Diedenhoven et al. (2014) geometric optics parameterization.

Single implementation of the parameterized single-scattering albedo (SSA)
and asymmetry parameter (g) calculations from:

    Van Diedenhoven, B., A. S. Ackerman, B. Cairns, and A. M. Fridlind (2014),
    A flexible parameterization for shortwave optical properties of ice crystals,
    J. Atmos. Sci., 71, 1763-1782.

Original supplementary code:
    https://www.researchgate.net/publication/259821840_ice_OP_parameterization

The parameterization uses geometric optics calculations (Macke et al., JAS,
1996) and applies to any convex particle given its aspect ratio, volume, and
projected area.  It is used by both the ice grain optical property generator
(geometric_optics_ice.py) and the bio-optical model (biooptical_funcs.py).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Parameterization tables (Figs. 4 and 7, Tables 2 and 3 of VD2014)
# ---------------------------------------------------------------------------

# SSA parameterization at ar=1
_a = [0.457593, 20.9738]

# SSA correction for AR != 1 (Table 2)
_NC1 = 3
_NC2 = 4
_c_ij = np.zeros((_NC1, _NC2, 2))
#   --- Plates ---
_c_ij[:, 0, 0] = [0.000527060, 0.309748, -2.58028]
_c_ij[:, 1, 0] = [0.00867596, -0.650188, -1.34949]
_c_ij[:, 2, 0] = [0.0382627, -0.198214, -0.674495]
_c_ij[:, 3, 0] = [0.0108558, -0.0356019, -0.141318]
#   --- Columns ---
_c_ij[:, 0, 1] = [0.000125752, 0.387729, -2.38400]
_c_ij[:, 1, 1] = [0.00797282, 0.456133, 1.29446]
_c_ij[:, 2, 1] = [0.00122800, -0.137621, -1.05868]
_c_ij[:, 3, 1] = [0.000212673, 0.0364655, 0.339646]

# Diffraction g parameterization
_b_gdiffr = [-0.822315, -1.20125, 0.996653]

# Raytracing g parameterization at ar=1
_p_a_eq_1 = [0.780550, 0.00510997, -0.0878268, 0.111549, -0.282453]

# g correction for AR != 1 (Table 3)
_NQ1 = 3
_NQ2 = 7
_q_ij = np.zeros((_NQ1, _NQ2, 2))
#   --- Plates ---
_q_ij[:, 0, 0] = [-0.00133106, -0.000782076, 0.00205422]
_q_ij[:, 1, 0] = [0.0408343, -0.00162734, 0.0240927]
_q_ij[:, 2, 0] = [0.525289, 0.418336, -0.818352]
_q_ij[:, 3, 0] = [0.443151, 1.53726, -2.40399]
_q_ij[:, 4, 0] = [0.00852515, 1.88625, -2.64651]
_q_ij[:, 5, 0] = [-0.123100, 0.983854, -1.29188]
_q_ij[:, 6, 0] = [-0.0376917, 0.187708, -0.235359]
#   --- Columns ---
_q_ij[:, 0, 1] = [-0.00189096, 0.000637430, 0.00157383]
_q_ij[:, 1, 1] = [0.00981029, 0.0409220, 0.00908004]
_q_ij[:, 2, 1] = [0.732647, 0.0539796, -0.665773]
_q_ij[:, 3, 1] = [-1.59927, -0.500870, 1.86375]
_q_ij[:, 4, 1] = [1.54047, 0.692547, -2.05390]
_q_ij[:, 5, 1] = [-0.707187, -0.374173, 1.01287]
_q_ij[:, 6, 1] = [0.125276, 0.0721572, -0.186466]

# Refractive index correction of asymmetry parameter
_c_g = np.zeros((2, 2))
_c_g[:, 0] = [0.96025050, 0.42918060]
_c_g[:, 1] = [0.94179149, -0.21600979]

# Correction for absorption
_s = [1.00014, 0.666094, -0.535922, -11.7454, 72.3600, -109.940]
_u = [-0.213038, 0.204016]

_DELTA = 0.3


def calc_ssa_and_g(ar, V, Area, reals, imags, wavelengths):
    """Calculate SSA and asymmetry parameter for a particle.

    Implements the Van Diedenhoven et al. (2014) parameterization.  The caller
    computes the geometry-dependent quantities (aspect ratio, volume, projected
    area) for their particle shape and passes them here.

    Args:
        ar: aspect ratio (depth/side_length for hex prisms,
            diameter/length for cylinders)
        V: particle volume (um^3)
        Area: mean projected area (um^2), typically total_surface_area / 4
        reals: real part of refractive index, array of length N
        imags: imaginary part of refractive index, array of length N
        wavelengths: wavelength array (um), array of length N

    Returns:
        ssa: single scattering albedo, array of length N
        g: asymmetry parameter, array of length N
    """
    col_pla = 1 if ar > 1.0 else 0
    log10_ar = np.log10(ar)

    ssa = np.empty(len(wavelengths))
    g = np.empty(len(wavelengths))

    for i in range(len(wavelengths)):
        mr = reals[i]
        mi = imags[i]
        wl = wavelengths[i]

        # --- Size parameters ---
        Chi_abs = mi / wl * V / Area  # absorption (Fig. 4, box 1)
        Chi_scat = 2.0 * np.pi * np.sqrt(Area / np.pi) / wl  # scattering (Fig. 7, box 1)

        # --- Single scattering albedo ---
        if Chi_abs > 0:
            # ar=1 (Fig. 4, box 2)
            w_1 = 1.0 - _a[0] * (1.0 - np.exp(-Chi_abs * _a[1]))
            l = np.zeros(_NC1)
            for j in range(_NC2):
                l[:] += _c_ij[:, j, col_pla] * log10_ar**j  # (Fig. 4, box 3)
            D_w = (
                l[0]
                * np.exp(-((np.log(Chi_abs) - l[2]) ** 2) / (2.0 * l[1] ** 2))
                / (Chi_abs * l[1] * np.sqrt(2.0 * np.pi))
            )  # (Fig. 4, box 3)
            w = w_1 + D_w  # (Fig. 4, box 4)
        else:
            w = 1.0

        # --- Asymmetry parameter ---
        # diffraction g (Fig. 7, box 2)
        g_diffr = (
            _b_gdiffr[0] * np.exp(_b_gdiffr[1] * np.log(Chi_scat)) + _b_gdiffr[2]
        )
        g_diffr = max(g_diffr, 0.5)

        # raytracing g at 862 nm (Fig. 7, box 3)
        g_1 = sum(_p_a_eq_1[j] * _DELTA**j for j in range(len(_p_a_eq_1)))

        # g correction for AR (Fig. 7, box 4)
        p_delta = np.zeros(_NQ1)
        for j in range(_NQ2):
            p_delta += _q_ij[:, j, col_pla] * log10_ar**j
        Dg = sum(p_delta[j] * _DELTA**j for j in range(_NQ1))
        g_rt = 2.0 * (g_1 + Dg) - 1.0  # (Fig. 7, box 5)

        # refractive index correction (Fig. 7, box 6)
        epsilon = _c_g[0, col_pla] + _c_g[1, col_pla] * log10_ar
        mr1 = 1.3038  # reference value @ 862 nm band
        # abs() added per corrigendum to original paper
        C_m = abs(
            (mr1 - epsilon) / (mr1 + epsilon) * (mr + epsilon) / (mr - epsilon)
        )

        # correction for absorption (Fig. 7, box 7)
        if Chi_abs > 0:
            C_w0 = sum(_s[j] * (1.0 - w) ** j for j in range(len(_s)))
            k_corr = log10_ar * _u[col_pla]
            C_w = C_w0 * (k_corr * w - k_corr + 1.0)
        else:
            C_w = 1.0

        # raytracing g at required wavelength (Fig. 7, box 9)
        g_rt_corr = g_rt * C_m * C_w

        # total asymmetry parameter, clamped to <= 1 (Fig. 7, box 9)
        g_tot = 1.0 / (2.0 * w) * ((2.0 * w - 1.0) * g_rt_corr + g_diffr)
        g_tot = min(g_tot, 1.0)

        ssa[i] = w
        g[i] = g_tot

    return ssa, g
