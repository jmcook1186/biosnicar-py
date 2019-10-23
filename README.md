# BioSNICAR_GO_PY

Translation of the BioSNICAR_GO model into Python (including translation of the original SNICAR model)

# Current Development status

21 Oct 2019: Mie scattering version of snicar functional. No unit testing done yet, pending access to previous SNICAR versions for benchmarking. Geometrical optics versions are not yet translated into Python and therefore setting GeometricOptics = 1 in the driver script raises an exception and will not run.

# In this repo

# Background

# How to use

# Known bugs and gotchas

1) Diffuse + Eddington
The combination of Eddington 2-stream approximation and diffuse incident radiation causes the albedo to go negative at wavelengths > ~1.4 um. Recommend using quadature or hemispheric mean approximations when diffuse
incident radiation is required.

2) SZA limits
The two stream aproximation seems to fall apart at high zenith angles (>~0.57). This is common to all versions of SNICAR and is explained in Toon et al. (1989).

3) Snow algae
While glacier algae MACs have been determined empirically, we have included only a hypothetical snow algae with potentially realistic pigemnt concentrations derived from the literature. A more accurate, empirically-derived set of single scattering optical properties for real snow algae is needed.


# Benchmarking

# Permissions

# Citation
