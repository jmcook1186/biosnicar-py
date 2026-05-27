Changelog
==========
All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Unreleased
----------

Fixed
~~~~~
- **Adding-doubling solver: incorrect albedo at SZA > ~55° for ice layers**
  (`#111 <https://github.com/jmcook1186/biosnicar-py/issues/111>`_).
  At oblique solar angles the direct beam can exceed the critical angle at
  ice's anomalous-dispersion bands (~2.93–3.09 µm), triggering total internal
  reflection (TIR) and producing a physically correct albedo of 1.0 in those
  bands. The Savitzky-Golay smoothing filter, applied as a post-processing
  step, treated these step discontinuities as noise and produced large ringing
  artefacts: values inside the TIR region were undershooting to 0.67–0.95, and
  the adjacent non-TIR transition bands were falsely elevated to 0.05–0.33.
  The fix identifies TIR bands before smoothing (``np.isclose(albedo, 1.0)``)
  and expands that mask outward by ``window_size // 2`` bands using
  ``scipy.ndimage.binary_dilation``. The raw solver output is preserved in
  this guard zone so the SG polynomial is never fitted across the
  discontinuity. The physical step at the TIR boundary—present in the Matlab
  reference (Whicker et al. 2022)—is reproduced correctly; the only
  transitions introduced at guard-zone edges are O(10⁻⁴) in magnitude and
  visually indistinguishable from the surrounding smoothed spectrum.

2.1.0 - (2023-07-20)
-------------

Changed
~~~~~~
- rename driver.py -> main.py
- remove read-the-docs in favour of custom site

Added
~~~~~~
- add wrapper func `get_albedo()` for one-line albedo prediction


