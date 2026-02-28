# Replace netCDF Optical Property Files with Compact Lookup Tables

## What changed

Replaced ~11,000 individual netCDF files (1.1 GB) with 11 compressed `.npz` lookup tables (129 MB), achieving an 88% storage reduction in the `Data/OP_data/480band/` directory while preserving bit-exact accuracy at all existing grid points.

## Files created

- `scripts/build_lookup_tables.py` — One-time conversion script (reads netCDF, writes `.npz`, validates bit-exact equality). Kept in repo so LUTs can be regenerated from any copy of the original netCDF data.
- `biosnicar/optical_properties/op_lookup.py` — `OpLookupTable` (1D radius index) and `HexOpLookupTable` (2D side/length index) classes with module-level caching and O(1) dict-based lookups.
- `Data/OP_data/480band/luts/*.npz` — 11 lookup table files:
  - `ice_sphere_Pic16.npz` (1,646 radii)
  - `ice_sphere_BH83_Pic16.npz` (1,508 radii)
  - `ice_sphere_BH83_Wrn08.npz` (1,471 radii)
  - `ice_sphere_BH83_Wrn84.npz` (1,471 radii)
  - `bubbly_air.npz` (549 radii)
  - `bubbly_water.npz` (549 radii)
  - `bubbly_air_BH83.npz` (1,535 radii)
  - `hex_Pic16.npz` (261 side/length pairs)
  - `hex_Wrn08.npz` (261 side/length pairs)
  - `hex_Wrn84.npz` (290 side/length pairs)
  - `water_sphere.npz` (1,646 radii)

## Files modified

- `biosnicar/classes/model_config.py` — Added `self.lut_dir` attribute pointing to `Data/OP_data/480band/luts/`.
- `biosnicar/optical_properties/column_OPs.py` — Replaced all `xr.open_dataset()` calls with LUT lookups. Removed `import xarray`. Added helper functions (`_sphere_lut_path`, `_hex_lut_path`, `_bubbly_air_lut_path`, `_bubbly_water_lut_path`, `_water_sphere_lut_path`) for LUT filename resolution from existing path config. Modified `add_water_coating()` to accept `ext_cff_mss_ice` as a parameter instead of reading it from a netCDF file.

## Files deleted

Six netCDF directories (~11,000 files, 1.1 GB total):
- `Data/OP_data/480band/ice_spherical_grains/`
- `Data/OP_data/480band/ice_spherical_grains_BH83/`
- `Data/OP_data/480band/bubbly_ice_files/`
- `Data/OP_data/480band/bubbly_ice_files_BH83/`
- `Data/OP_data/480band/ice_hexagonal_columns/`
- `Data/OP_data/480band/water_spherical_grains/`

## Verification

- All 11 LUTs validated bit-exact against source netCDF files (0 mismatches).
- Full test suite: 242/242 passed, including AD solver benchmark (3,000 sims, 1e-5 tolerance vs Matlab), config fuzzer (252 combos), and var fuzzer (108 combos).
- Tests passed both before and after deleting the netCDF directories.
- `main.py` smoke test produces valid output (BBA = 0.4186).

## No API changes

YAML configuration, class constructors, and `get_albedo.get()` all work identically. The LUT filename is derived automatically from the existing `sphere_ice_path`, `bubbly_ice_path`, etc. config strings, so tests that override paths to BH83 variants continue to work.

---

# Phase 2: Migrate remaining netCDF files (LAP, FSDS, rfidx, Fresnel)

## What changed

Converted all remaining netCDF files to `.npz` format: 40 LAP (light absorbing particle) files, 647 FSDS (illumination) files, and 2 standalone files (`rfidx_ice.nc`, `fl_reflection_diffuse.nc`). Total: 689 `.nc` files removed, replaced by 689 `.npz` files.

## Storage reduction

| Directory | Before (.nc) | After (.npz) |
|---|---|---|
| `lap/` | ~3.8 MB | 596 KB |
| `fsds/` | ~26 MB | 2.7 MB |
| `rfidx_ice` | 66 KB | 32 KB |
| `fl_reflection_diffuse` | 45 KB | 24 KB |
| **Total 480band/** | **~1.13 GB** | **133 MB** |

Overall reduction: **~88%** across both phases.

## Files created

- `scripts/build_lap_fsds_npz.py` — Conversion script for LAP, FSDS, rfidx_ice, and fl_reflection_diffuse files. Validates bit-exact equality (with NaN-aware comparison) against source files.
- `Data/OP_data/480band/lap/*.npz` — 40 files (variables: `ss_alb`, `asm_prm`, plus `ext_cff_mss`/`ext_xsc`/`ext_cff_mss_ncl` where present in the original).
- `Data/OP_data/480band/fsds/*.npz` — 647 files (all original data_vars preserved).
- `Data/OP_data/480band/rfidx_ice.npz` — Refractive index data (8 data_vars + `wvl` coordinate).
- `Data/OP_data/480band/fl_reflection_diffuse.npz` — Fresnel reflectance data (6 coordinate variables).

## Files modified

- `biosnicar/classes/impurity.py` — Replaced `xr.open_dataset()` with `np.load()`. File extension `.nc` → `.npz` via `os.path.splitext`. Removed `import xarray`.
- `biosnicar/classes/illumination.py` — Replaced `xr.open_dataset()` with `np.load()`. Changed `.nc` extension to `.npz` in path construction. Used `.copy()` on `flx_frc_sfc` (mutated in-place). Removed `import xarray`.
- `biosnicar/classes/ice.py` — Replaced `xr.open_dataset()` with `np.load()` for both `rfidx_ice` and `fl_reflection_diffuse`. Used `.copy()` on `ref_idx_im` (mutated in-place by CDOM code). Removed `import xarray`.
- `biosnicar/optical_properties/mie_coated_water_spheres.py` — Replaced `xr.open_dataset()` with `np.load()`. Removed `import xarray`.
- `biosnicar/optical_properties/geometric_optics_ice.py` — Updated path from `rfidx_ice.nc` to `rfidx_ice.npz`, replaced xarray with numpy load.
- `biosnicar/optical_properties/ssps_spheres_generator.py` — Updated `rfidx_ice.nc` to `rfidx_ice.npz`, replaced xarray with numpy load.
- `biosnicar/inputs.yaml` — Changed `FN_ICE` path from `.nc` to `.npz`.
- `app_inputs.yaml` — Changed `FN_ICE` path from `.nc` to `.npz`.

## Files deleted

- 40 `.nc` files in `Data/OP_data/480band/lap/`
- 647 `.nc` files in `Data/OP_data/480band/fsds/`
- `Data/OP_data/480band/rfidx_ice.nc`
- `Data/OP_data/480band/fl_reflection_diffuse.nc`

## Key implementation details

- **npz read-only arrays**: `.npz` arrays are memory-mapped and read-only. `.copy()` is required for arrays that are mutated in-place (`ref_idx_im` for CDOM, `flx_slr` for zero-replacement).
- **NaN handling**: `rfidx_ice` contains NaN values in `im_co2ice`, and 11 dust LAP files have NaN in `ext_cff_mss_ncl`. Validation uses `np.array_equal(equal_nan=True)`.
- **wvl coordinate in rfidx_ice.npz**: `mie_coated_water_spheres.py` reads `temp["wvl"]` from rfidx_ice, so the `wvl` dimension coordinate is explicitly included in the npz file.

## Verification

- All 689 files validated bit-exact against source netCDF files (0 mismatches after NaN-aware comparison fix).
- Full test suite: 242/242 passed after both code changes and file deletion.
- `main.py` smoke test produces valid output.
