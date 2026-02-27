#!/usr/bin/env python3
"""Convert LAP, FSDS, and standalone netCDF files to .npz format.

Converts:
  - Data/OP_data/480band/lap/*.nc → Data/OP_data/480band/lap/*.npz
  - Data/OP_data/480band/fsds/*.nc → Data/OP_data/480band/fsds/*.npz
  - Data/OP_data/480band/rfidx_ice.nc → Data/OP_data/480band/rfidx_ice.npz
  - Data/OP_data/480band/fl_reflection_diffuse.nc → Data/OP_data/480band/fl_reflection_diffuse.npz

Usage:
    python scripts/build_lap_fsds_npz.py
"""

import sys
import numpy as np
import xarray as xr
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data" / "OP_data" / "480band"


def convert_lap_files():
    """Convert each LAP netCDF to a minimal .npz with only runtime variables.

    Each LAP file is used by the Impurity class which reads:
      - ss_alb (always)
      - asm_prm (always)
      - ext_cff_mss (for standard impurities)
      - ext_xsc (for algae: name == 'ga' or 'sa')
      - ext_cff_mss_ncl (for coated impurities)

    We store all of these that exist in each file so any loading mode works.
    """
    lap_dir = DATA_DIR / "lap"
    nc_files = sorted(lap_dir.glob("*.nc"))
    print(f"\n[LAP] Converting {len(nc_files)} files in {lap_dir}")

    results = []
    for nc_path in nc_files:
        with xr.open_dataset(str(nc_path)) as ds:
            arrays = {}
            arrays["ss_alb"] = ds["ss_alb"].values.copy()
            arrays["asm_prm"] = ds["asm_prm"].values.copy()

            # Store all mac variants that exist in this file
            for var in ["ext_cff_mss", "ext_xsc", "ext_cff_mss_ncl"]:
                if var in ds:
                    arrays[var] = ds[var].values.copy()

        npz_path = nc_path.with_suffix(".npz")
        np.savez_compressed(str(npz_path), **arrays)
        results.append((nc_path, npz_path, arrays))
        print(f"  {nc_path.name} -> {npz_path.name} ({len(arrays)} vars)")

    return results


def validate_lap_files(results):
    """Validate LAP .npz files against source netCDF."""
    errors = 0
    for nc_path, npz_path, _ in results:
        lut = np.load(str(npz_path))
        with xr.open_dataset(str(nc_path)) as ds:
            for var in lut.files:
                orig = ds[var].values
                stored = lut[var]
                if not np.array_equal(orig, stored, equal_nan=True):
                    print(f"  MISMATCH: {nc_path.name} var={var}")
                    errors += 1
    return errors


def convert_fsds_files():
    """Convert each FSDS netCDF to an individual .npz.

    Each FSDS file contains flx_frc_sfc (480,). We do a simple 1:1 conversion
    for each file, keeping the same name with .npz extension.
    """
    fsds_dir = DATA_DIR / "fsds"
    nc_files = sorted(fsds_dir.glob("*.nc"))
    print(f"\n[FSDS] Converting {len(nc_files)} files in {fsds_dir}")

    results = []
    for nc_path in nc_files:
        with xr.open_dataset(str(nc_path)) as ds:
            arrays = {}
            for var in ds.data_vars:
                arrays[var] = ds[var].values.copy()

        npz_path = nc_path.with_suffix(".npz")
        np.savez_compressed(str(npz_path), **arrays)
        results.append((nc_path, npz_path, arrays))

    print(f"  Converted {len(results)} files")
    return results


def validate_fsds_files(results):
    """Validate FSDS .npz files against source netCDF."""
    errors = 0
    for nc_path, npz_path, _ in results:
        lut = np.load(str(npz_path))
        with xr.open_dataset(str(nc_path)) as ds:
            for var in lut.files:
                orig = ds[var].values
                stored = lut[var]
                if not np.array_equal(orig, stored, equal_nan=True):
                    print(f"  MISMATCH: {nc_path.name} var={var}")
                    errors += 1
    return errors


def convert_rfidx_ice():
    """Convert rfidx_ice.nc to rfidx_ice.npz.

    Variables: im_Pic16, re_Pic16, im_Wrn84, re_Wrn84, im_Wrn08, re_Wrn08,
               im_co2ice, re_co2ice
    """
    nc_path = DATA_DIR / "rfidx_ice.nc"
    npz_path = DATA_DIR / "rfidx_ice.npz"
    print(f"\n[rfidx_ice] Converting {nc_path}")

    with xr.open_dataset(str(nc_path)) as ds:
        arrays = {}
        # Include the wvl coordinate (used by mie_coated_water_spheres)
        arrays["wvl"] = ds["wvl"].values.copy()
        for var in ds.data_vars:
            arrays[var] = ds[var].values.copy()

    np.savez_compressed(str(npz_path), **arrays)
    print(f"  -> {npz_path.name} ({len(arrays)} vars: {list(arrays.keys())})")
    return nc_path, npz_path, arrays


def validate_rfidx_ice(nc_path, npz_path, arrays):
    """Validate rfidx_ice.npz against source."""
    errors = 0
    lut = np.load(str(npz_path))
    with xr.open_dataset(str(nc_path)) as ds:
        for var in lut.files:
            if var == "wvl":
                orig = ds["wvl"].values
            else:
                orig = ds[var].values
            stored = lut[var]
            if not np.array_equal(orig, stored, equal_nan=True):
                print(f"  MISMATCH: rfidx_ice var={var}")
                errors += 1
    return errors


def convert_fl_reflection():
    """Convert fl_reflection_diffuse.nc to fl_reflection_diffuse.npz.

    This file stores Fresnel reflectance data as coordinates (not data variables).
    Variables: R_dif_fa_ice_Wrn84, R_dif_fb_ice_Wrn84, R_dif_fa_ice_Wrn08,
               R_dif_fb_ice_Wrn08, R_dif_fa_ice_Pic16, R_dif_fb_ice_Pic16
    """
    nc_path = DATA_DIR / "fl_reflection_diffuse.nc"
    npz_path = DATA_DIR / "fl_reflection_diffuse.npz"
    print(f"\n[fl_reflection_diffuse] Converting {nc_path}")

    with xr.open_dataset(str(nc_path)) as ds:
        arrays = {}
        # These are stored as coordinates in the netCDF, not data variables
        for coord_name in ds.coords:
            if coord_name == "wvl":
                continue  # skip the dimension coordinate
            arrays[coord_name] = ds.coords[coord_name].values.copy()

    np.savez_compressed(str(npz_path), **arrays)
    print(f"  -> {npz_path.name} ({len(arrays)} vars: {list(arrays.keys())})")
    return nc_path, npz_path, arrays


def validate_fl_reflection(nc_path, npz_path, arrays):
    """Validate fl_reflection_diffuse.npz against source."""
    errors = 0
    lut = np.load(str(npz_path))
    with xr.open_dataset(str(nc_path)) as ds:
        for var in lut.files:
            orig = ds.coords[var].values
            stored = lut[var]
            if not np.array_equal(orig, stored, equal_nan=True):
                print(f"  MISMATCH: fl_reflection_diffuse var={var}")
                errors += 1
    return errors


def main():
    total_errors = 0

    print("=" * 60)
    print("Converting LAP, FSDS, and standalone netCDF files to .npz")
    print("=" * 60)

    # LAP files
    lap_results = convert_lap_files()

    # FSDS files
    fsds_results = convert_fsds_files()

    # Standalone files
    rfidx_nc, rfidx_npz, rfidx_arrays = convert_rfidx_ice()
    fl_nc, fl_npz, fl_arrays = convert_fl_reflection()

    # Validation
    print("\n" + "=" * 60)
    print("Validating all conversions")
    print("=" * 60)

    print(f"\n  Validating {len(lap_results)} LAP files...")
    errors = validate_lap_files(lap_results)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    print(f"\n  Validating {len(fsds_results)} FSDS files...")
    errors = validate_fsds_files(fsds_results)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    print(f"\n  Validating rfidx_ice.npz...")
    errors = validate_rfidx_ice(rfidx_nc, rfidx_npz, rfidx_arrays)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    print(f"\n  Validating fl_reflection_diffuse.npz...")
    errors = validate_fl_reflection(fl_nc, fl_npz, fl_arrays)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    # Summary
    print("\n" + "=" * 60)
    npz_count = (
        len(lap_results) + len(fsds_results) + 2  # rfidx + fl_reflection
    )
    print(f"Converted {npz_count} files total")
    print(f"Validation errors: {total_errors}")

    if total_errors > 0:
        print("\nFAILED: Some values did not match source files!")
        sys.exit(1)
    else:
        print("\nSUCCESS: All values are bit-exact with source netCDF files.")
        sys.exit(0)


if __name__ == "__main__":
    main()
