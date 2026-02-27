#!/usr/bin/env python3
"""Build compact .npz lookup tables from netCDF optical property files.

Run this script ONCE while the original netCDF files still exist.
It reads every netCDF file per category, extracts only the variables
used at runtime, consolidates them into compressed .npz archives,
and validates bit-exact equality against the source files.

Usage:
    python scripts/build_lookup_tables.py
"""

import os
import re
import sys
import numpy as np
import xarray as xr
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data" / "OP_data" / "480band"
LUT_DIR = DATA_DIR / "luts"


def extract_radius_from_filename(fname, prefix):
    """Extract integer radius from filename like 'ice_Pic16_0200.nc' or 'bbl_0200.nc'."""
    stem = Path(fname).stem  # e.g. 'ice_Pic16_0200' or 'bbl_0200'
    suffix = stem[len(prefix):]
    return int(suffix)


def build_sphere_lut(src_dir, ri_stub, lut_name):
    """Build LUT for spherical ice grains.

    Variables extracted: ss_alb, ext_cff_mss, asm_prm, ext_cff_vlm, sca_cff_vlm

    Args:
        src_dir: directory containing netCDF files
        ri_stub: e.g. 'ice_Pic16_' - filename prefix before the radius
        lut_name: output filename (without path)
    """
    print(f"  Building {lut_name} from {src_dir}")
    files = sorted(Path(src_dir).glob(f"{ri_stub}*.nc"))
    if not files:
        print(f"    WARNING: No files found matching {ri_stub}*.nc in {src_dir}")
        return

    radii = []
    ss_alb_list = []
    ext_cff_mss_list = []
    asm_prm_list = []
    ext_cff_vlm_list = []
    sca_cff_vlm_list = []

    for f in files:
        r = extract_radius_from_filename(f.name, ri_stub)
        radii.append(r)
        with xr.open_dataset(str(f)) as ds:
            ss_alb_list.append(ds["ss_alb"].values.copy())
            ext_cff_mss_list.append(ds["ext_cff_mss"].values.copy())
            asm_prm_list.append(ds["asm_prm"].values.copy())
            ext_cff_vlm_list.append(ds["ext_cff_vlm"].values.copy())
            sca_cff_vlm_list.append(ds["sca_cff_vlm"].values.copy())

    radii = np.array(radii, dtype=np.int32)
    sort_idx = np.argsort(radii)
    radii = radii[sort_idx]

    out_path = LUT_DIR / lut_name
    np.savez_compressed(
        str(out_path),
        radii=radii,
        ss_alb=np.array(ss_alb_list)[sort_idx],
        ext_cff_mss=np.array(ext_cff_mss_list)[sort_idx],
        asm_prm=np.array(asm_prm_list)[sort_idx],
        ext_cff_vlm=np.array(ext_cff_vlm_list)[sort_idx],
        sca_cff_vlm=np.array(sca_cff_vlm_list)[sort_idx],
    )
    print(f"    Wrote {len(radii)} entries to {out_path}")
    return out_path, files, ri_stub


def build_bubbly_air_lut(src_dir, lut_name):
    """Build LUT for bubbly ice (air bubbles).

    Variables extracted: sca_cff_vlm, asm_prm
    """
    print(f"  Building {lut_name} from {src_dir}")
    files = sorted(Path(src_dir).glob("bbl_[0-9]*.nc"))
    if not files:
        print(f"    WARNING: No bbl_*.nc files found in {src_dir}")
        return

    radii = []
    sca_cff_vlm_list = []
    asm_prm_list = []

    for f in files:
        r = extract_radius_from_filename(f.name, "bbl_")
        radii.append(r)
        with xr.open_dataset(str(f)) as ds:
            sca_cff_vlm_list.append(ds["sca_cff_vlm"].values.copy())
            asm_prm_list.append(ds["asm_prm"].values.copy())

    radii = np.array(radii, dtype=np.int32)
    sort_idx = np.argsort(radii)
    radii = radii[sort_idx]

    out_path = LUT_DIR / lut_name
    np.savez_compressed(
        str(out_path),
        radii=radii,
        sca_cff_vlm=np.array(sca_cff_vlm_list)[sort_idx],
        asm_prm=np.array(asm_prm_list)[sort_idx],
    )
    print(f"    Wrote {len(radii)} entries to {out_path}")
    return out_path, files


def build_bubbly_water_lut(src_dir, lut_name):
    """Build LUT for bubbly ice (water inclusions).

    Variables extracted: sca_cff_vlm, ext_cff_vlm, asm_prm
    """
    print(f"  Building {lut_name} from {src_dir}")
    files = sorted(Path(src_dir).glob("bbl_water_*.nc"))
    if not files:
        print(f"    WARNING: No bbl_water_*.nc files found in {src_dir}")
        return

    radii = []
    sca_cff_vlm_list = []
    ext_cff_vlm_list = []
    asm_prm_list = []

    for f in files:
        r = extract_radius_from_filename(f.name, "bbl_water_")
        radii.append(r)
        with xr.open_dataset(str(f)) as ds:
            sca_cff_vlm_list.append(ds["sca_cff_vlm"].values.copy())
            ext_cff_vlm_list.append(ds["ext_cff_vlm"].values.copy())
            asm_prm_list.append(ds["asm_prm"].values.copy())

    radii = np.array(radii, dtype=np.int32)
    sort_idx = np.argsort(radii)
    radii = radii[sort_idx]

    out_path = LUT_DIR / lut_name
    np.savez_compressed(
        str(out_path),
        radii=radii,
        sca_cff_vlm=np.array(sca_cff_vlm_list)[sort_idx],
        ext_cff_vlm=np.array(ext_cff_vlm_list)[sort_idx],
        asm_prm=np.array(asm_prm_list)[sort_idx],
    )
    print(f"    Wrote {len(radii)} entries to {out_path}")
    return out_path, files


def build_hex_lut(src_dir, ri_stub, lut_name):
    """Build LUT for hexagonal ice columns (2D index: side x length).

    Variables extracted: ss_alb, ext_cff_mss, asm_prm
    """
    print(f"  Building {lut_name} from {src_dir}")
    files = sorted(Path(src_dir).glob(f"{ri_stub}*.nc"))
    if not files:
        print(f"    WARNING: No files found matching {ri_stub}*.nc in {src_dir}")
        return

    sides = []
    lengths = []
    ss_alb_list = []
    ext_cff_mss_list = []
    asm_prm_list = []

    for f in files:
        stem = f.stem  # e.g. 'ice_Pic16_10000_10000'
        parts = stem[len(ri_stub):].split("_")
        side = int(parts[0])
        length = int(parts[1])
        sides.append(side)
        lengths.append(length)
        with xr.open_dataset(str(f)) as ds:
            ss_alb_list.append(ds["ss_alb"].values.copy())
            ext_cff_mss_list.append(ds["ext_cff_mss"].values.copy())
            asm_prm_list.append(ds["asm_prm"].values.copy())

    sides = np.array(sides, dtype=np.int32)
    lengths = np.array(lengths, dtype=np.int32)
    # Sort by (side, length) for consistent ordering
    sort_idx = np.lexsort((lengths, sides))
    sides = sides[sort_idx]
    lengths = lengths[sort_idx]

    out_path = LUT_DIR / lut_name
    np.savez_compressed(
        str(out_path),
        sides=sides,
        lengths=lengths,
        ss_alb=np.array(ss_alb_list)[sort_idx],
        ext_cff_mss=np.array(ext_cff_mss_list)[sort_idx],
        asm_prm=np.array(asm_prm_list)[sort_idx],
    )
    print(f"    Wrote {len(sides)} entries to {out_path}")
    return out_path, files, ri_stub


def build_water_sphere_lut(src_dir, lut_name):
    """Build LUT for water spherical grains.

    Variables extracted: asm_prm, ext_cff_vlm, sca_cff_vlm
    """
    print(f"  Building {lut_name} from {src_dir}")
    files = sorted(Path(src_dir).glob("water_grain_*.nc"))
    if not files:
        print(f"    WARNING: No water_grain_*.nc files found in {src_dir}")
        return

    radii = []
    asm_prm_list = []
    ext_cff_vlm_list = []
    sca_cff_vlm_list = []

    for f in files:
        r = extract_radius_from_filename(f.name, "water_grain_")
        radii.append(r)
        with xr.open_dataset(str(f)) as ds:
            asm_prm_list.append(ds["asm_prm"].values.copy())
            ext_cff_vlm_list.append(ds["ext_cff_vlm"].values.copy())
            sca_cff_vlm_list.append(ds["sca_cff_vlm"].values.copy())

    radii = np.array(radii, dtype=np.int32)
    sort_idx = np.argsort(radii)
    radii = radii[sort_idx]

    out_path = LUT_DIR / lut_name
    np.savez_compressed(
        str(out_path),
        radii=radii,
        asm_prm=np.array(asm_prm_list)[sort_idx],
        ext_cff_vlm=np.array(ext_cff_vlm_list)[sort_idx],
        sca_cff_vlm=np.array(sca_cff_vlm_list)[sort_idx],
    )
    print(f"    Wrote {len(radii)} entries to {out_path}")
    return out_path, files


def validate_sphere_lut(lut_path, files, prefix):
    """Validate sphere LUT against source netCDF files."""
    lut = np.load(str(lut_path))
    radii = lut["radii"]
    radius_to_idx = {int(r): i for i, r in enumerate(radii)}
    errors = 0

    for f in files:
        r = extract_radius_from_filename(f.name, prefix)
        idx = radius_to_idx[r]
        with xr.open_dataset(str(f)) as ds:
            for var in ["ss_alb", "ext_cff_mss", "asm_prm", "ext_cff_vlm", "sca_cff_vlm"]:
                orig = ds[var].values
                stored = lut[var][idx]
                if not np.array_equal(orig, stored):
                    print(f"    MISMATCH: {f.name} var={var}")
                    errors += 1
    return errors


def validate_bubbly_air_lut(lut_path, files):
    """Validate bubbly air LUT against source netCDF files."""
    lut = np.load(str(lut_path))
    radii = lut["radii"]
    radius_to_idx = {int(r): i for i, r in enumerate(radii)}
    errors = 0

    for f in files:
        r = extract_radius_from_filename(f.name, "bbl_")
        idx = radius_to_idx[r]
        with xr.open_dataset(str(f)) as ds:
            for var in ["sca_cff_vlm", "asm_prm"]:
                orig = ds[var].values
                stored = lut[var][idx]
                if not np.array_equal(orig, stored):
                    print(f"    MISMATCH: {f.name} var={var}")
                    errors += 1
    return errors


def validate_bubbly_water_lut(lut_path, files):
    """Validate bubbly water LUT against source netCDF files."""
    lut = np.load(str(lut_path))
    radii = lut["radii"]
    radius_to_idx = {int(r): i for i, r in enumerate(radii)}
    errors = 0

    for f in files:
        r = extract_radius_from_filename(f.name, "bbl_water_")
        idx = radius_to_idx[r]
        with xr.open_dataset(str(f)) as ds:
            for var in ["sca_cff_vlm", "ext_cff_vlm", "asm_prm"]:
                orig = ds[var].values
                stored = lut[var][idx]
                if not np.array_equal(orig, stored):
                    print(f"    MISMATCH: {f.name} var={var}")
                    errors += 1
    return errors


def validate_hex_lut(lut_path, files, prefix):
    """Validate hex column LUT against source netCDF files."""
    lut = np.load(str(lut_path))
    sides = lut["sides"]
    lengths = lut["lengths"]
    key_to_idx = {(int(s), int(l)): i for i, (s, l) in enumerate(zip(sides, lengths))}
    errors = 0

    for f in files:
        stem = f.stem
        parts = stem[len(prefix):].split("_")
        side = int(parts[0])
        length = int(parts[1])
        idx = key_to_idx[(side, length)]
        with xr.open_dataset(str(f)) as ds:
            for var in ["ss_alb", "ext_cff_mss", "asm_prm"]:
                orig = ds[var].values
                stored = lut[var][idx]
                if not np.array_equal(orig, stored):
                    print(f"    MISMATCH: {f.name} var={var}")
                    errors += 1
    return errors


def validate_water_sphere_lut(lut_path, files):
    """Validate water sphere LUT against source netCDF files."""
    lut = np.load(str(lut_path))
    radii = lut["radii"]
    radius_to_idx = {int(r): i for i, r in enumerate(radii)}
    errors = 0

    for f in files:
        r = extract_radius_from_filename(f.name, "water_grain_")
        idx = radius_to_idx[r]
        with xr.open_dataset(str(f)) as ds:
            for var in ["asm_prm", "ext_cff_vlm", "sca_cff_vlm"]:
                orig = ds[var].values
                stored = lut[var][idx]
                if not np.array_equal(orig, stored):
                    print(f"    MISMATCH: {f.name} var={var}")
                    errors += 1
    return errors


def main():
    LUT_DIR.mkdir(parents=True, exist_ok=True)
    total_errors = 0
    build_results = []

    print("=" * 60)
    print("Building compact lookup tables from netCDF files")
    print("=" * 60)

    # 1. Ice spheres - default (Pic16 only)
    print("\n[1/11] Ice spheres (default, Pic16)")
    result = build_sphere_lut(
        DATA_DIR / "ice_spherical_grains" / "ice_Pic16",
        "ice_Pic16_",
        "ice_sphere_Pic16.npz",
    )
    if result:
        build_results.append(("sphere", result))

    # 2-4. Ice spheres - BH83 (Pic16, Wrn08, Wrn84)
    for i, ri in enumerate(["Pic16", "Wrn08", "Wrn84"], start=2):
        print(f"\n[{i}/11] Ice spheres (BH83, {ri})")
        result = build_sphere_lut(
            DATA_DIR / "ice_spherical_grains_BH83" / f"ice_{ri}",
            f"ice_{ri}_",
            f"ice_sphere_BH83_{ri}.npz",
        )
        if result:
            build_results.append(("sphere", result))

    # 5. Bubbly air (default)
    print("\n[5/11] Bubbly air (default)")
    result = build_bubbly_air_lut(
        DATA_DIR / "bubbly_ice_files",
        "bubbly_air.npz",
    )
    if result:
        build_results.append(("bubbly_air", result))

    # 6. Bubbly water (default)
    print("\n[6/11] Bubbly water (default)")
    result = build_bubbly_water_lut(
        DATA_DIR / "bubbly_ice_files",
        "bubbly_water.npz",
    )
    if result:
        build_results.append(("bubbly_water", result))

    # 7. Bubbly air (BH83)
    print("\n[7/11] Bubbly air (BH83)")
    result = build_bubbly_air_lut(
        DATA_DIR / "bubbly_ice_files_BH83",
        "bubbly_air_BH83.npz",
    )
    if result:
        build_results.append(("bubbly_air", result))

    # 8-10. Hex columns (Pic16, Wrn08, Wrn84)
    for i, ri in enumerate(["Pic16", "Wrn08", "Wrn84"], start=8):
        print(f"\n[{i}/11] Hex columns ({ri})")
        result = build_hex_lut(
            DATA_DIR / "ice_hexagonal_columns" / f"ice_{ri}",
            f"ice_{ri}_",
            f"hex_{ri}.npz",
        )
        if result:
            build_results.append(("hex", result))

    # 11. Water spheres
    print("\n[11/11] Water spheres")
    result = build_water_sphere_lut(
        DATA_DIR / "water_spherical_grains",
        "water_sphere.npz",
    )
    if result:
        build_results.append(("water_sphere", result))

    # Validation
    print("\n" + "=" * 60)
    print("Validating all LUTs against source netCDF files")
    print("=" * 60)

    for kind, result in build_results:
        if kind == "sphere":
            lut_path, files, prefix = result
            print(f"\n  Validating {lut_path.name} ({len(files)} files)...")
            errors = validate_sphere_lut(lut_path, files, prefix)
        elif kind == "bubbly_air":
            lut_path, files = result
            print(f"\n  Validating {lut_path.name} ({len(files)} files)...")
            errors = validate_bubbly_air_lut(lut_path, files)
        elif kind == "bubbly_water":
            lut_path, files = result
            print(f"\n  Validating {lut_path.name} ({len(files)} files)...")
            errors = validate_bubbly_water_lut(lut_path, files)
        elif kind == "hex":
            lut_path, files, prefix = result
            print(f"\n  Validating {lut_path.name} ({len(files)} files)...")
            errors = validate_hex_lut(lut_path, files, prefix)
        elif kind == "water_sphere":
            lut_path, files = result
            print(f"\n  Validating {lut_path.name} ({len(files)} files)...")
            errors = validate_water_sphere_lut(lut_path, files)
        total_errors += errors
        if errors == 0:
            print(f"    OK - all values match")

    # Summary
    print("\n" + "=" * 60)
    lut_files = list(LUT_DIR.glob("*.npz"))
    total_size = sum(f.stat().st_size for f in lut_files)
    print(f"Created {len(lut_files)} LUT files in {LUT_DIR}")
    print(f"Total size: {total_size / 1024 / 1024:.1f} MB")
    print(f"Validation errors: {total_errors}")

    if total_errors > 0:
        print("\nFAILED: Some values did not match source files!")
        sys.exit(1)
    else:
        print("\nSUCCESS: All LUT values are bit-exact with source netCDF files.")
        sys.exit(0)


if __name__ == "__main__":
    main()
