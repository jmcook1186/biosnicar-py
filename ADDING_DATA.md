# Adding New Data to BioSNICAR

This guide covers how to add new data to the three main data categories:
light-absorbing particles (LAP), incoming irradiance (FSDS), and ice grain
lookup tables (LUT).

All data lives under `Data/OP_data/480band/`. Each category uses a single
compressed `.npz` archive loaded once at runtime.

---

## 1. Light-Absorbing Particles (LAP)

**Archive:** `Data/OP_data/480band/lap.npz`
**Build script:** `scripts/build_consolidated_npz.py`

### Key scheme

Each impurity contributes several arrays keyed as `{stem}__{variable}`:

| Variable          | Shape    | Description                    |
|-------------------|----------|--------------------------------|
| `ss_alb`          | (480,)   | Single-scattering albedo       |
| `asm_prm`         | (480,)   | Asymmetry parameter            |
| `ext_cff_mss`     | (480,)   | Mass extinction coefficient    |
| `ext_xsc`         | (480,)   | Extinction cross-section (algae) |
| `ext_cff_mss_ncl` | (480,)   | Coated mass extinction coeff.  |

Not all variables exist for every impurity. A `_names` key stores all
impurity stems.

### Steps to add a new impurity

1. Generate optical properties using Mie theory or geometric optics (see
   `biosnicar/optical_properties/`). The output should be a `.npz` file with
   at least `ss_alb`, `asm_prm`, and one of the `ext_*` variables above,
   each of shape `(480,)`.

2. Place the individual `.npz` file in a temporary `Data/OP_data/480band/lap/`
   directory alongside the other individual files (you can extract them from
   `lap.npz` if needed, or just add yours).

3. Run the consolidation script:
   ```
   python scripts/build_consolidated_npz.py
   ```
   This rebuilds `lap.npz` and validates bit-exact equality.

4. Add a config entry in `biosnicar/inputs.yaml` under `IMPURITIES`:
   ```yaml
   MY_IMP:
     NAME: "my_imp"
     FILE: "my_impurity_file.npz"
     COATED: False
     UNIT: 0        # 0 = ppb, 1 = cells/mL
     CONC: [0, 0]
   ```

5. If your impurity uses `ext_xsc` (cell-based), set `UNIT: 1` and use
   `NAME: "sa"` or `NAME: "ga"` (the name controls which mac variable is
   read by `Impurity.__init__`).

6. If your impurity has coated properties (`ext_cff_mss_ncl`), set
   `COATED: True`.

---

## 2. Incoming Irradiance (FSDS)

**Archive:** `Data/OP_data/480band/fsds.npz`
**Build script:** `scripts/build_consolidated_npz.py`

### Key scheme

Each illumination file is stored as `{filename_stem}` -> `ndarray(480,)`.
The variable stored is `flx_frc_sfc` (fractional surface flux). A `_names`
key stores all stems.

Keys follow the pattern:
`swnb_480bnd_{atmosphere}_{cloud}` (diffuse) or
`swnb_480bnd_{atmosphere}_clr_SZA{angle}` (direct).

### Steps to add a new illumination spectrum

1. Generate an FSDS file (e.g., from a radiative transfer model). It should
   contain a `flx_frc_sfc` variable of shape `(480,)` in a `.npz` file.

2. Place the individual `.npz` in a temporary `Data/OP_data/480band/fsds/`
   directory.

3. Run the consolidation script:
   ```
   python scripts/build_consolidated_npz.py
   ```

4. Add a new stub to `ILLUMINATION_FILE_STUBS` in `biosnicar/inputs.yaml`:
   ```yaml
   ILLUMINATION_FILE_STUBS: ["swnb_480bnd_mlw", ..., "swnb_480bnd_mynew"]
   ```
   The stub index maps to the `INCOMING` integer selector.

---

## 3. Ice Grain Lookup Tables (LUT)

**Archive directory:** `Data/OP_data/480band/luts/`
**Build script:** `scripts/build_lookup_tables.py`

The LUT system stores optical properties for ice grains and bubbles indexed
by effective radius. Each refractive-index variant gets its own `.npz` file.

### LUT files

| File                            | Contents                              |
|---------------------------------|---------------------------------------|
| `sphere_ice_{ri}.npz`          | Spherical ice grains (per RI variant) |
| `bubbly_ice_{ri}.npz`          | Bubbly ice (per RI variant)           |
| `hex_ice.npz`                  | Hexagonal ice columns                 |
| `sphere_water.npz`             | Water-coated spheres                  |

### Steps to add new ice grain optical properties

1. Generate per-radius netCDF files using the Mie/GO code in
   `biosnicar/optical_properties/`. Each file should contain `ss_alb`,
   `asm_prm`, and `ext_cff_mss` (all shape `(480,)`).

2. Place them in the appropriate source directory (e.g.,
   `Data/OP_data/480band/ice_spherical_grains/`).

3. Run the LUT build script:
   ```
   python scripts/build_lookup_tables.py
   ```
   This consolidates individual files into the compact `.npz` archives in
   `Data/OP_data/480band/luts/` and validates all values.

---

## General Notes

- All arrays must have exactly 480 elements (one per wavelength band).
- The build scripts validate bit-exact equality against source files.
- After rebuilding any archive, run the full test suite: `pytest tests/`
- The `_names` keys in LAP and FSDS archives list all available entries
  for discoverability.
