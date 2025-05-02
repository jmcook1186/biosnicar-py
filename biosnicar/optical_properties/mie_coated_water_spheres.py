#!/usr/bin/python

"""
Copyright (C) 2020  Niklas Bohn (GFZ, <nbohn@gfz-potsdam.de>),
German Research Centre for Geosciences (GFZ, <https://www.gfz-potsdam.de>)
"""
import numpy as np
import pandas as pd
import xarray as xr
import miepython as mie
from scipy.interpolate import interp1d
from scipy.special import jv, yv
from tqdm import tqdm


def fill_nans_scipy1(padata, pkind="nearest"):
    """Interpolates data to fill nan values.

    Args:
    padata: source data with np.NaN values
    pkind:  kind of interpolation (see scipy.interpolate.interp1d docs)

    Returns:
        f(aindexes): data with interpolated values instead of nans
    """

    aindexes = np.arange(padata.shape[0])
    (agood_indexes,) = np.where(np.isfinite(padata))
    f = interp1d(
        agood_indexes,
        padata[agood_indexes],
        bounds_error=False,
        copy=False,
        fill_value="extrapolate",
        kind=pkind,
    )

    return f(aindexes)


def miecoated_ab3(m1, m2, x, y):
    """Computation of Mie Coefficients.

    Computes a_n, b_n, of orders n=1 to nmax, complex
    refractive index m=m'+im", and size parameters
    x=k0*a, y=k0*b where k0 = wave number in the ambient medium for
    coated spheres, a = inner radius, b = outer radius
    m1, m2 = inner, outer refractive index;

    p. 183 in Bohren and Huffman (1983) BEWI:TDD122 but
    using the bottom equation on p. 483 for chi_prime (Matzler 2002).

    Args:
        m1: refractive index for inner sphere
        m2: refractive index for outher sphere
        x:  size parameter for inner sphere
        y:  size parameter for outer sphere

    Returns:
       Mie Coefficients a_n and b_n
    """

    M = m2 / m1
    nmax = int(round(2 + y + 4 * y ** (1 / 3)))
    n = np.arange(1, nmax + 1, 1)
    nu = n + 0.5
    u = m1 * x
    v = m2 * x
    w = m2 * y
    su = np.sqrt(0.5 * np.pi * u)
    sv = np.sqrt(0.5 * np.pi * v)
    sw = np.sqrt(0.5 * np.pi * w)
    sy = np.sqrt(0.5 * np.pi * y)
    pu = su * jv(nu, u)
    py = sy * jv(nu, y)
    pv = sv * jv(nu, v)
    pw = sw * jv(nu, w)
    p1u = np.concatenate([np.array([np.sin(u)]), pu[0 : nmax - 1]])
    p1y = np.concatenate([np.array([np.sin(y)]), py[0 : nmax - 1]])
    p1v = np.concatenate([np.array([np.sin(v)]), pv[0 : nmax - 1]])
    p1w = np.concatenate([np.array([np.sin(w)]), pw[0 : nmax - 1]])
    ppw = p1w - n * pw / w
    ppy = p1y - n * py / y
    chv = -sv * yv(nu, v)
    chw = -sw * yv(nu, w)
    chy = -sy * yv(nu, y)
    ch1y = np.concatenate([np.array([np.cos(y)]), chy[0 : nmax - 1]])
    gsy = py - 1j * chy
    gs1y = p1y - 1j * ch1y
    gspy = gs1y - n * gsy / y
    du = p1u / pu - n / u
    dv = p1v / pv - n / v
    dw = p1w / pw - n / w
    chpw = chw * dw - 1.0 / pw
    uu = M * du - dv
    vv = du / M - dv
    pvi = 1.0 / pv
    aaa = pv * uu / (chv * uu + pvi)
    bbb = pv * vv / (chv * vv + pvi)
    aa1 = ppw - aaa * chpw
    aa2 = pw - aaa * chw
    bb1 = ppw - bbb * chpw
    bb2 = pw - bbb * chw
    aa = (py * aa1 - m2 * ppy * aa2) / (gsy * aa1 - m2 * gspy * aa2)
    bb = (m2 * py * bb1 - ppy * bb2) / (m2 * gsy * bb1 - gspy * bb2)

    f = np.concatenate((aa, bb))
    f = f.reshape(2, nmax)

    return f


def miecoated(m1, m2, x, y):
    """Mie Efficiencies of coated spheres.

    Calculates Mie efficiencies for given complex refractive-index
    ratios m1=m1'+im1", m2= m2'+im2" of kernel
    and coating, resp., and size parameters x=k0*a, y=k0*b where
    k0 = wave number in ambient  medium, a,b = inner,
    outer sphere radius, using complex Mie Coefficients an and bn for
    n=1 to nmax, s. Bohren and Huffman (1983)
    BEWI:TDD122, p. 181-185,483.

    opt selects the function "Miecoated_ab.." for an and bn, n=1 to nmax.

    Note that 0<=x<=y (Matzler, 2002).

    Args:
        m1:  refractive-index ratio of kernel
        m2:  refractive-index ratio of coating
        x:   size parameter for inner sphere
        y:   size parameter for outer sphere

    Returns:
        qext: extinction efficiency
        qsca: scatterign efficiency
        qb: backscattering efficiency
        asy: asymmetry parameter
        qratio: qb/qsca

    """

    # to reduce computing time
    if x == y:
        return mie(m1, y)

    # to avoid a singularity at x=0
    elif x == 0:
        return mie(m2, y)

    # to reduce computing time
    elif m1 == m2:
        return mie(m1, y)

    # this is the normal situation
    elif x > 0:
        nmax = int(round(2 + y + 4 * y ** (1 / 3)))
        n1 = nmax - 1
        n = np.arange(1, nmax + 1, 1)
        cn = 2 * n + 1
        c1n = n * (n + 2) / (n + 1)
        c2n = cn / n / (n + 1)
        y2 = y * y

        f = miecoated_ab3(m1, m2, x, y)

        anp = f[0, :].real
        anpp = f[0, :].imag
        bnp = f[1, :].real
        bnpp = f[1, :].imag

        # displaced numbers used for asymmetry parameter, p. 120
        g1 = np.zeros((4, nmax))
        g1[0, :n1] = anp[1 : nmax + 1]
        g1[1, :n1] = anpp[1 : nmax + 1]
        g1[2, :n1] = bnp[1 : nmax + 1]
        g1[3, :n1] = bnpp[1 : nmax + 1]

        dn = cn * (anp + bnp)
        q = sum(dn)
        qext = 2 * q / y2
        en = cn * (anp * anp + anpp * anpp + bnp * bnp + bnpp * bnpp)
        q = sum(en)
        qsca = 2 * q / y2
        qabs = qext - qsca
        fn = (f[0, :] - f[1, :]) * cn
        gn = (-1) ^ n
        f_3 = fn * gn
        q = sum(f_3)
        qb = q * q.conj() / y2
        asy1 = c1n * (
            anp * g1[0, :] + anpp * g1[1, :] + bnp * g1[2, :] + bnpp * g1[3, :]
        )
        asy2 = c2n * (anp * bnp + anpp * bnpp)
        asy = 4 / y2 * sum(asy1 + asy2) / qsca
        qratio = qb / qsca

        return qext, qsca, qabs, qb, asy, qratio


def miecoated_driver(rice, rwater, fn_ice, rf_ice, fn_water, wvl):
    """Driver for miecoated.

    Originally written by Christian Matzler
    (see Matzler, 2002). The driver convolves the efficiency factors with the particle
    dimensions to return the cross sections for extinction, scattering and
    absorption plus the asymmetry parameter, q ratio and single scattering albedo.

    Note that the code includes an interpolation regime. This is because the original
    code produced NaNs for a few wavelengths at certain size parameters, particularly
    in the mid NIR wavelengths.

    Adapted from Matlab code by Joseph Cook, University of Sheffield, UK (2017).

    Args:
        rice:     inner sphere diameter in microns
        rwater:   outer sphere diameter in microns (i.e. total coated sphere,
                        not water layer thickness)
        fn_ice:   path to csv file containing refractive index of ice (Warren, 1984)
        fn_water: path to csv file containing refractive index of liquid water
                    (Segelstein, 1981)
        wvl:      wavelength which should be calculated (in microns)

        Returns:
            res: tuple containing cross sections for extinction, scattering and
            absorption plus the asymmetry parameter, q ratio and single scattering
            albedo
    """

    # calculate volume and density of sphere

    # cross sectional area of ice core
    XSArea_inner = np.pi * (rice**2)

    # cross sectional area of water layer sphere
    XSArea_outer = np.pi * (rwater**2) - XSArea_inner
    TotalXS = np.pi * (rwater**2)

    # density of water at 1 degree C in kg m-3
    WatDensity = 999

    # density of ice in kg m-3
    IceDensity = 934

    IceVol = 4 / 3 * np.pi * rice**3
    WatVol = (4 / 3 * np.pi * rwater**3) - IceVol
    TotalVol = 4 / 3 * np.pi * rwater**3

    IceMass = IceVol * IceDensity
    WatMass = WatVol * WatDensity
    TotalMass = IceMass + WatMass

    # read in refractive indices of ice and liquid water

    temp = xr.open_dataset(fn_ice)
    ref_index_ice = np.zeros(shape=(2, len(wvl)))
    wvl_ice = temp["wvl"].values

    if rf_ice == 0:
        n_ice = temp["re_Wrn84"].values
        k_ice = temp["im_Wrn84"].values

    if rf_ice == 1:
        n_ice = temp["re_Wrn08"].values
        k_ice = temp["im_Wrn08"].values

    if rf_ice == 2:
        n_ice = temp["re_Pic16"].values
        k_ice = temp["im_Pic16"].values

    ref_index_water = pd.read_csv(fn_water)
    wvl_water = np.zeros(ref_index_water.shape[0])
    n_water = np.zeros(ref_index_water.shape[0])
    k_water = np.zeros(ref_index_water.shape[0])

    for ii in range(ref_index_water.shape[0]):
        wvl_water[ii] = ref_index_water.at[ii, "wl"]
        n_water[ii] = ref_index_water.at[ii, "n"]
        k_water[ii] = ref_index_water.at[ii, "k"]

    n_water_interp = np.interp(x=wvl, xp=wvl_water, fp=n_water)
    k_water_interp = np.interp(x=wvl, xp=wvl_water, fp=k_water)

    extinction = np.zeros(len(wvl))
    scattering = np.zeros(len(wvl))
    absorption = np.zeros(len(wvl))
    backscattering = np.zeros(len(wvl))
    asymmetry = np.zeros(len(wvl))
    q_ratio = np.zeros(len(wvl))
    ssa = np.zeros(len(wvl))

    for ii in tqdm(range(len(wvl))):
        # size parameters for inner and outer spheres
        x = 2 * np.pi * rice / (wvl[ii])
        y = 2 * np.pi * rwater / (wvl[ii])

        m1 = n_ice[ii] + k_ice[ii] * 1j
        m2 = n_water_interp[ii] + k_water_interp[ii] * 1j

        # call Miecoated and return efficiencies
        qext, qsca, qabs, qb, asy, qratio = miecoated(m1=m1, m2=m2, x=x, y=y)

        # append efficiencies to lists
        extinction[ii] = qext
        scattering[ii] = qsca
        absorption[ii] = qabs
        backscattering[ii] = qb
        asymmetry[ii] = asy
        q_ratio[ii] = qratio
        ssa[ii] = qsca / qext

    # replace any possible nans with values estimated by cubic interpolation
    extinction = fill_nans_scipy1(padata=extinction, pkind="nearest")
    scattering = fill_nans_scipy1(padata=scattering, pkind="nearest")
    absorption = fill_nans_scipy1(padata=absorption, pkind="nearest")
    backscattering = fill_nans_scipy1(padata=backscattering, pkind="nearest")
    asymmetry = fill_nans_scipy1(padata=asymmetry, pkind="nearest")
    q_ratio = fill_nans_scipy1(padata=q_ratio, pkind="nearest")
    ssa = fill_nans_scipy1(padata=ssa, pkind="nearest")

    # calculate cross sections from efficiency factors
    ExtXC = extinction * TotalXS
    ScaXC = scattering * TotalXS
    AbsXC = absorption * TotalXS

    ExtXCvol = extinction * TotalVol
    ScaXCvol = scattering * TotalVol
    AbsXCvol = absorption * TotalVol

    ExtXCmass = extinction * TotalMass
    ScaXCmass = scattering * TotalMass
    AbsXCmass = absorption * TotalMass

    # print fraction of sphere made up of water by mass and volume
    water_frac_mss = 100 - (IceMass / TotalMass) * 100
    water_frac_vol = 100 - (IceVol / TotalVol) * 100
    ice_frac_mss = 100 - water_frac_mss
    ice_frac_vol = 100 - water_frac_vol

    # density of particle as average weighted by mass of components
    part_dens = (IceDensity * ice_frac_mss / 100) + (WatDensity * water_frac_mss / 100)

    res = {
        "extinction": extinction,
        "scattering": scattering,
        "absorption": absorption,
        "backscattering": backscattering,
        "asymmetry": asymmetry,
        "q_ratio": q_ratio,
        "ssa": ssa,
        "extinction_cross_section": ExtXC,
        "scattering_cross_section": ScaXC,
        "absorption_cross_section": AbsXC,
        "extinction_volume_cross_section": ExtXCvol,
        "scattering_volume_cross_section": ScaXCvol,
        "absorption_volume_cross_section": AbsXCvol,
        "extinction_mass_cross_section": ExtXCmass,
        "scattering_mass_cross_section": ScaXCmass,
        "absorption_mass_cross_section": AbsXCmass,
        "water_mass_fraction": water_frac_mss,
        "water_volume_fraction": water_frac_vol,
        "ice_mass_fraction": ice_frac_mss,
        "ice_volume_fraction": ice_frac_vol,
        "particle_density": part_dens,
    }

    return res


if __name__ == "__main__":
    pass
