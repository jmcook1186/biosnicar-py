#!/usr/bin/env python3
"""Demonstrate the parameter sweep API.

Runs a solar zenith angle x grain radius sweep and prints summary results.

Usage:
    python scripts/sweep_demo.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from biosnicar.drivers.sweep import parameter_sweep


def main():
    # --- Example 1: SZA x grain radius sweep ---
    print("Running SZA x grain radius sweep (20 combinations)...")
    df = parameter_sweep(
        params={
            "solzen": [30, 40, 50, 60, 70],
            "rds": [100, 200, 500, 1000],
        },
        progress=True,
    )

    print("\nBroadband albedo pivot table (SZA vs grain radius):\n")
    pivot = df.pivot_table(values="BBA", index="solzen", columns="rds")
    print(pivot.to_string(float_format="{:.4f}".format))

    # --- Example 2: Impurity concentration sweep ---
    print("\n\nRunning black carbon concentration sweep...")
    df_bc = parameter_sweep(
        params={"impurity.0.conc": [0, 100, 1000, 10000, 100000]},
        progress=True,
    )

    print("\nBBA vs black carbon concentration:\n")
    for _, row in df_bc.iterrows():
        print(f"  BC = {row['impurity.0.conc']:>8.0f} ppb  ->  BBA = {row['BBA']:.4f}")

    # --- Example 3: Spectral output ---
    print("\n\nRunning single point with spectral output...")
    df_spec = parameter_sweep(
        params={"solzen": [50]},
        return_spectral=True,
        progress=False,
    )

    albedo = df_spec.iloc[0]["albedo"]
    print(f"  Spectral albedo: {len(albedo)} bands, "
          f"range [{albedo.min():.4f}, {albedo.max():.4f}]")


if __name__ == "__main__":
    main()
