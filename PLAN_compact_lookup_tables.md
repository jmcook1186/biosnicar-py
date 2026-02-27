# Plan: Replace netCDF Optical Property Files with Compact Lookup Tables

## Context

biosnicar-py stores ~12,000 pre-computed Mie theory results as individual netCDF files (1.1 GB total in `Data/OP_data/480band/`). Each file is indexed by grain radius and contains 480-wavelength optical property spectra. At runtime, `column_OPs.py:get_layer_OPs()` opens one file per layer per simulation, extracts 2-5 arrays of shape `(480,)`, and uses them in the radiative transfer solver. The netCDF format adds ~95% overhead (metadata, unused variables, per-file size distribution grids). The actual data consumed at runtime is ~120 MB of raw float64 arrays.

**Goal**: Replace ~12,000 netCDF files with compact `.npz` lookup tables (~50-70 MB compressed), achieving a **94-96% storage reduction** while preserving **bit-exact accuracy** at all existing grid points. The netCDF optical property files are fully removed — no fallback, no conditional logic.

## Recommended Approach: Compact Lookup Tables (sole code path)

**Why this approach over alternatives:**
- **On-the-fly Mie calculation**: Rejected — the files contain lognormal-averaged properties from 5,000+ Mie evaluations per radius (see `ssps_spheres_generator.py`); hex columns use geometric optics, not Mie. Would be slow and cannot guarantee bit-exact match.
- **SVD/PCA compression**: Rejected — lossy; truncation error is hard to bound and could affect scientific outputs.
- **Neural network**: Rejected — introduces approximation error incompatible with benchmark tests requiring 1e-5 tolerance.

## Files to Create

### 1. `scripts/build_lookup_tables.py` — One-time conversion script
Reads every netCDF file per category, extracts only the variables used at runtime, and consolidates into compressed `.npz` archives. Includes built-in validation that verifies bit-exact equality against source files after writing. This script is run ONCE before any code changes, while the netCDF files still exist.

**Output files** (stored in `Data/OP_data/480band/luts/`):

| LUT file | Source | Index | Variables |
|---|---|---|---|
| `ice_sphere_Pic16.npz` | `ice_spherical_grains/ice_Pic16/` | radii (1D) | ss_alb, ext_cff_mss, asm_prm, ext_cff_vlm, sca_cff_vlm |
| `ice_sphere_BH83_Pic16.npz` | `ice_spherical_grains_BH83/ice_Pic16/` | radii | same |
| `ice_sphere_BH83_Wrn08.npz` | `ice_spherical_grains_BH83/ice_Wrn08/` | radii | same |
| `ice_sphere_BH83_Wrn84.npz` | `ice_spherical_grains_BH83/ice_Wrn84/` | radii | same |
| `bubbly_air.npz` | `bubbly_ice_files/bbl_*.nc` | radii | sca_cff_vlm, asm_prm |
| `bubbly_water.npz` | `bubbly_ice_files/bbl_water_*.nc` | radii | sca_cff_vlm, ext_cff_vlm, asm_prm |
| `bubbly_air_BH83.npz` | `bubbly_ice_files_BH83/bbl_*.nc` | radii | sca_cff_vlm, asm_prm |
| `hex_Pic16.npz` | `ice_hexagonal_columns/ice_Pic16/` | sides + lengths (2D) | ss_alb, ext_cff_mss, asm_prm |
| `hex_Wrn08.npz` | `ice_hexagonal_columns/ice_Wrn08/` | sides + lengths | same |
| `hex_Wrn84.npz` | `ice_hexagonal_columns/ice_Wrn84/` | sides + lengths | same |
| `water_sphere.npz` | `water_spherical_grains/` | radii | asm_prm, ext_cff_vlm, sca_cff_vlm |

### 2. `biosnicar/optical_properties/op_lookup.py` — LUT loader class
- `OpLookupTable`: loads an `.npz`, builds `{radius: index}` dict for O(1) lookup, provides `get(radius, var_name) -> ndarray(480,)`.
- `HexOpLookupTable`: subclass for 2D (side, length) indexing.
- Class-level `_cache` dict prevents redundant disk reads across layers/simulations.

## Files to Modify

### 3. `biosnicar/optical_properties/column_OPs.py` — Core change
**Replace** all `xr.open_dataset()` calls with LUT lookups (no fallback, no conditionals):

- **layer_type 0, shp < 4** (lines 61-96): Replace netCDF open with ice sphere LUT `get()`.
- **layer_type 0, shp == 4** (lines 51-59): Replace with hex column LUT `get(side, length, var)`.
- **layer_type 1-2** (lines 100-178): Replace with bubbly ice LUT `get()`.
- **layer_type 3** (lines 181-227): Replace with ice sphere + water sphere LUT `get()`.
- **add_water_coating** (line 276): Replace netCDF open with LUT `get()` for `ext_cff_mss`.
- Remove `import xarray as xr` from this file (no longer needed here).

The LUT filename is derived from `model_config.sphere_ice_path` / `bubbly_ice_path` etc., so when tests override these paths to BH83 variants, the correct BH83 LUT is resolved automatically.

### 4. `biosnicar/classes/model_config.py` — Add LUT directory path
Add `self.lut_dir` attribute pointing to `Data/OP_data/480band/luts/`. No `USE_LUT` toggle needed since the LUT is now the only code path.

### 5. `biosnicar/inputs.yaml` and `app_inputs.yaml` — Update PATHS
- Replace `SPHERE_ICE`, `SPHERE_WATER`, `HEX_ICE`, `BUBBLY_ICE` path entries with references to the LUT directory.
- Or: keep the path config entries for LUT filename resolution (the path string is used to derive which `.npz` to load, e.g. "BH83" in the path → `*_BH83_*.npz`).

### 6. Delete netCDF optical property directories
After conversion and validation, remove from the repository:
- `Data/OP_data/480band/ice_spherical_grains/` (1,646 files)
- `Data/OP_data/480band/ice_spherical_grains_BH83/` (4,450 files)
- `Data/OP_data/480band/bubbly_ice_files/` (1,098 files)
- `Data/OP_data/480band/bubbly_ice_files_BH83/` (1,535 files)
- `Data/OP_data/480band/ice_hexagonal_columns/` (812 files)
- `Data/OP_data/480band/water_spherical_grains/` (1,646 files)

**Files kept as-is** (small, different purpose):
- `Data/OP_data/480band/lap/` — 40 LAP files, ~1 MB total (loaded by `Impurity` class)
- `Data/OP_data/480band/fsds/` — illumination files (loaded by `Illumination` class)
- `Data/OP_data/480band/rfidx_ice.nc` — 65 KB refractive index (loaded by `Ice` class)
- `Data/OP_data/480band/fl_reflection_diffuse.nc` — 44 KB Fresnel data (loaded by `Ice` class)
- Various CSV reference files (CDOM, water RI, surface reflectance, wavelengths)

## Complete File Selection Dimensions Audit

Every user-configurable parameter that determines which optical property file gets loaded has been verified:

### Parameters that SELECT which LUT + index to use:

1. **LAYER_TYPE** (per layer): 0=granular→sphere/hex LUT, 1-2=solid→bubbly LUT, 3=mixed→sphere+water LUTs
2. **SHP** (grain shape, layer_type 0 only): 0-3→sphere LUT (with analytical asphericity correction for 1-3), 4→hex column LUT
3. **RF** (refractive index): 0=Wrn84, 1=Wrn08, 2=Pic16 → selects which `.npz` file. Note: default solver only has Pic16 spheres; BH83 has all three.
4. **Solver variant** (default vs BH83): configured via SPHERE_ICE and BUBBLY_ICE paths → determines LUT filename prefix.
5. **RDS** (radius): integer μm, used as lookup key within the LUT
6. **HEX_SIDE + HEX_LENGTH** (hex prisms only): 2D lookup key in hex column LUT

### Parameters that modify data AFTER loading (no LUT selection impact):
- SHP_FCTR, AR → asphericity correction applied analytically to loaded g
- WATER_COATING → triggers runtime Mie calc; only ext_cff_mss still loaded from LUT
- RHO, DZ, LWC, LWC_PCT_BBL → used in mixing/optical depth calculations
- CDOM → loads tiny CSV, not part of the LUT system

### All 11 LUT files cover every existing directory:

| LUT | Solver | RI | Category |
|---|---|---|---|
| `ice_sphere_Pic16.npz` | Default | Pic16 | Only RI in default dir |
| `ice_sphere_BH83_{Pic16,Wrn08,Wrn84}.npz` | BH83 | All 3 | All 3 subdirs |
| `bubbly_air.npz` | Default | N/A | RI not in file path |
| `bubbly_water.npz` | Default | N/A | Only in default dir |
| `bubbly_air_BH83.npz` | BH83 | N/A | Air only (no BH83 water files exist) |
| `hex_{Pic16,Wrn08,Wrn84}.npz` | N/A | All 3 | All 3 subdirs |
| `water_sphere.npz` | N/A | N/A | No RI or solver variant |

## Verification

1. **Conversion validation**: `build_lookup_tables.py` verifies every value in every `.npz` against the source netCDF after writing (run before deleting netCDF files).
2. **Existing tests must pass unchanged**: The benchmark tests (`test_AD_solver`, `test_compare_pyBBA_to_matBBA`, etc.) use BH83 data and require 1e-5 tolerance vs Matlab — these must pass with LUTs.
3. **Fuzzer tests**: Config fuzzer (252 combos) and var fuzzer (108 combos) must all pass.
4. **Quick smoke test**: Run `main.py` to verify default config produces valid output.

## Implementation Order

1. `scripts/build_lookup_tables.py` — write and run to generate `.npz` files (while netCDFs still exist)
2. Validate all `.npz` files against source netCDFs (built into the script)
3. `biosnicar/optical_properties/op_lookup.py` — loader class
4. Modify `biosnicar/classes/model_config.py` — add `lut_dir` attribute
5. Modify `biosnicar/optical_properties/column_OPs.py` — replace all `xr.open_dataset` with LUT lookups
6. Run full existing test suite to confirm no regressions
7. Delete the 6 netCDF optical property directories (~11,187 files, 1.1 GB)
8. Run tests once more to confirm everything works without the netCDF files

## Risk Mitigation

- **Validation before deletion**: netCDF files are only deleted after both the conversion script AND the full test suite confirm correctness.
- **No API changes**: `get_albedo.get()`, `setup_snicar()`, class constructors — all unchanged.
- **Conversion script preserved**: `scripts/build_lookup_tables.py` remains in the repo so LUTs can be regenerated if needed from any copy of the original netCDF data.
- **Git LFS**: The `.npz` files (~50-70 MB total) should be tracked with Git LFS.

## User Interface

No changes. All YAML configuration parameters work identically. The only difference is internal: numpy array lookups instead of netCDF file opens.
