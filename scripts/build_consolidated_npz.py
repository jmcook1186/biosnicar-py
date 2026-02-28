#!/usr/bin/env python3
"""Consolidate individual LAP and FSDS .npz files into single archives.

Produces:
  - Data/OP_data/480band/lap.npz   (40 impurity files → one archive)
  - Data/OP_data/480band/fsds.npz  (647 illumination files → one archive)

Key schemes:
  LAP:  {stem}__{variable}  e.g. "bc_ChCB_rn40_dns1270__ss_alb"
        Plus a "_names" key listing all impurity stems.
  FSDS: {stem}              e.g. "swnb_480bnd_hmn_clr_SZA50"
        The only variable stored is flx_frc_sfc (the sole runtime variable).
        Plus a "_names" key listing all filename stems.

Usage:
    python scripts/build_consolidated_npz.py
"""

import sys
import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data" / "OP_data" / "480band"


def build_lap():
    """Consolidate individual LAP .npz files into a single lap.npz."""
    lap_dir = DATA_DIR / "lap"
    files = sorted(lap_dir.glob("*.npz"))
    print(f"\n[LAP] Found {len(files)} files in {lap_dir}")

    arrays = {}
    names = []

    for f in files:
        stem = f.stem
        names.append(stem)
        data = np.load(str(f))
        for var in data.files:
            key = f"{stem}__{var}"
            arrays[key] = data[var]

    arrays["_names"] = np.array(names)

    out_path = DATA_DIR / "lap.npz"
    np.savez_compressed(str(out_path), **arrays)
    print(f"  -> {out_path.name} ({len(arrays)} keys, {len(names)} impurities)")
    return out_path, files


def build_fsds():
    """Consolidate individual FSDS .npz files into a single fsds.npz."""
    fsds_dir = DATA_DIR / "fsds"
    files = sorted(fsds_dir.glob("*.npz"))
    print(f"\n[FSDS] Found {len(files)} files in {fsds_dir}")

    arrays = {}
    names = []

    for f in files:
        stem = f.stem
        names.append(stem)
        data = np.load(str(f))
        # Most files have flx_frc_sfc; the TOA file has flx_frc_toa
        if "flx_frc_sfc" in data.files:
            arrays[stem] = data["flx_frc_sfc"]
        elif "flx_frc_toa" in data.files:
            arrays[stem] = data["flx_frc_toa"]
        else:
            raise ValueError(f"No expected flux key in {f.name}: {data.files}")

    arrays["_names"] = np.array(names)

    out_path = DATA_DIR / "fsds.npz"
    np.savez_compressed(str(out_path), **arrays)
    print(f"  -> {out_path.name} ({len(arrays)} keys, {len(names)} illumination files)")
    return out_path, files


def validate_lap(out_path, source_files):
    """Validate every value in lap.npz against individual source files."""
    consolidated = np.load(str(out_path))
    errors = 0

    for f in source_files:
        stem = f.stem
        original = np.load(str(f))
        for var in original.files:
            key = f"{stem}__{var}"
            if key not in consolidated:
                print(f"  MISSING KEY: {key}")
                errors += 1
                continue
            if not np.array_equal(original[var], consolidated[key], equal_nan=True):
                print(f"  MISMATCH: {key}")
                errors += 1

    return errors


def validate_fsds(out_path, source_files):
    """Validate every value in fsds.npz against individual source files."""
    consolidated = np.load(str(out_path))
    errors = 0

    for f in source_files:
        stem = f.stem
        original = np.load(str(f))
        if "flx_frc_sfc" in original.files:
            orig_val = original["flx_frc_sfc"]
        else:
            orig_val = original["flx_frc_toa"]
        if stem not in consolidated:
            print(f"  MISSING KEY: {stem}")
            errors += 1
            continue
        if not np.array_equal(orig_val, consolidated[stem], equal_nan=True):
            print(f"  MISMATCH: {stem}")
            errors += 1

    return errors


def main():
    total_errors = 0

    print("=" * 60)
    print("Consolidating LAP and FSDS into single .npz archives")
    print("=" * 60)

    lap_path, lap_files = build_lap()
    fsds_path, fsds_files = build_fsds()

    print("\n" + "=" * 60)
    print("Validating consolidated archives")
    print("=" * 60)

    print(f"\n  Validating lap.npz against {len(lap_files)} source files...")
    errors = validate_lap(lap_path, lap_files)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    print(f"\n  Validating fsds.npz against {len(fsds_files)} source files...")
    errors = validate_fsds(fsds_path, fsds_files)
    total_errors += errors
    print(f"    {'OK' if errors == 0 else f'ERRORS: {errors}'}")

    print("\n" + "=" * 60)
    if total_errors > 0:
        print(f"FAILED: {total_errors} validation errors!")
        sys.exit(1)
    else:
        print("SUCCESS: All values are bit-exact with source files.")
        sys.exit(0)


if __name__ == "__main__":
    main()
