# BioSNICAR_GO_PY

Translation of the BioSNICAR_GO model into Python (including translation of the original SNICAR model)

# In this repo

# Background

# How to use

# Known bugs and gotchas

1) This version is currently limited to 5 layer configurations. This is because of the way the A, B, D and E are solved near line 425 - 460. In the FORTRAN version and in the Toon et al. paper, there is a simple if-else condition that
applies a variation of the calculations depending upon whether the layer number is odd or even. Converting to Python means switching from a 1-indexing scheme to a 0-indexing scheme which breaks that pattern and causes the equations to call indexes that exceed the max no layers. In this release, a workaround is to semi-hard code the layer numbers into the if-else statements (i.e. rather than "if i is even, then..., we use if i = 2 or 4, then...).

2) SZA limits
The two stream aproximation seems to fall apart at high zenith angles (>~0.57). This is common to all versions of SNICAR and is explained in Toon et al. (1989).

3)Snow algae
While glacir algae MACs have been determined empirically, we have included only a hypothetical snow algae with potentially
realistic pigemnt concentrations derived from the literature. 
A more accurate, empirically-derived set of single scattering
optical properties for real snow algae is needed.


# Benchmarking

# Permissions

# Citation
