#!/usr/bin/env python3
"""Demonstrate the parameter sweep API.

Examples cover 2-parameter, 3-parameter, and impurity sweeps with optional
matplotlib plotting.

Usage:
    python scripts/sweep_demo.py            # text output only
    python scripts/sweep_demo.py --plot     # also produce figures
"""

import argparse

import numpy as np

from biosnicar.drivers.sweep import parameter_sweep


def main():
    parser = argparse.ArgumentParser(description="Parameter sweep demo")
    parser.add_argument(
        "--plot", action="store_true", help="Show matplotlib figures"
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Example 1: Two-parameter sweep (SZA x grain radius)
    # ------------------------------------------------------------------
    print("=== Example 1: SZA x grain radius (20 combinations) ===\n")
    df = parameter_sweep(
        params={
            "solzen": [30, 40, 50, 60, 70],
            "rds": [100, 200, 500, 1000],
        },
        progress=True,
    )

    pivot = df.pivot_table(values="BBA", index="solzen", columns="rds")
    print("\nBroadband albedo pivot table (SZA vs grain radius):\n")
    print(pivot.to_string(float_format="{:.4f}".format))

    if args.plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        for rds in pivot.columns:
            ax.plot(pivot.index, pivot[rds], marker="o", label=f"r = {rds} µm")
        ax.set_xlabel("Solar zenith angle (°)")
        ax.set_ylabel("Broadband albedo")
        ax.set_title("BBA vs SZA for different grain radii")
        ax.legend()
        fig.tight_layout()

    # ------------------------------------------------------------------
    # Example 2: Impurity concentration sweep
    # ------------------------------------------------------------------
    print("\n\n=== Example 2: Black carbon concentration sweep ===\n")
    df_bc = parameter_sweep(
        params={"impurity.0.conc": [0, 100, 1000, 10000, 100000]},
        progress=True,
    )

    print("\nBBA vs black carbon concentration:\n")
    for _, row in df_bc.iterrows():
        print(f"  BC = {row['impurity.0.conc']:>8.0f} ppb  ->  BBA = {row['BBA']:.4f}")

    if args.plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.semilogx(
            df_bc["impurity.0.conc"].replace(0, 1),  # avoid log(0)
            df_bc["BBA"],
            marker="s",
            color="k",
        )
        ax.set_xlabel("Black carbon concentration (ppb)")
        ax.set_ylabel("Broadband albedo")
        ax.set_title("Albedo darkening from black carbon")
        fig.tight_layout()

    # ------------------------------------------------------------------
    # Example 3: Three-parameter sweep (SZA x density x grain radius)
    # ------------------------------------------------------------------
    print("\n\n=== Example 3: SZA x density x grain radius (36 combinations) ===\n")
    df3 = parameter_sweep(
        params={
            "solzen": [40, 50, 60],
            "rho": [400, 600, 800],
            "rds": [200, 500, 1000, 5000],
        },
        progress=True,
    )

    # Summarise: for each SZA, show the rho x rds BBA table
    for sza in sorted(df3["solzen"].unique()):
        subset = df3[df3["solzen"] == sza]
        piv = subset.pivot_table(values="BBA", index="rho", columns="rds")
        print(f"\n  SZA = {sza}°  —  BBA (density vs grain radius):\n")
        print("  " + piv.to_string(float_format="{:.4f}".format).replace("\n", "\n  "))

    if args.plot:
        import matplotlib.pyplot as plt

        # Heatmap for each SZA
        szas = sorted(df3["solzen"].unique())
        fig, axes = plt.subplots(1, len(szas), figsize=(5 * len(szas), 4), sharey=True)
        if len(szas) == 1:
            axes = [axes]

        for ax, sza in zip(axes, szas):
            subset = df3[df3["solzen"] == sza]
            piv = subset.pivot_table(values="BBA", index="rho", columns="rds")
            im = ax.imshow(
                piv.values,
                aspect="auto",
                origin="lower",
                vmin=df3["BBA"].min(),
                vmax=df3["BBA"].max(),
            )
            ax.set_xticks(range(len(piv.columns)))
            ax.set_xticklabels(piv.columns)
            ax.set_yticks(range(len(piv.index)))
            ax.set_yticklabels(piv.index)
            ax.set_xlabel("Grain radius (µm)")
            ax.set_title(f"SZA = {sza}°")

        axes[0].set_ylabel("Density (kg/m³)")
        fig.colorbar(im, ax=axes, label="Broadband albedo", shrink=0.8)
        fig.suptitle("BBA across SZA, density, and grain radius", y=1.02)
        fig.tight_layout()

    # ------------------------------------------------------------------
    # Example 4: Spectral output with plotting
    # ------------------------------------------------------------------
    print("\n\n=== Example 4: Spectral albedo for three grain radii ===\n")
    df_spec = parameter_sweep(
        params={"rds": [200, 1000, 5000]},
        return_spectral=True,
        progress=False,
    )

    wavelengths = np.arange(0.205, 4.999, 0.01)
    for _, row in df_spec.iterrows():
        alb = row["albedo"]
        print(
            f"  rds = {row['rds']:>5.0f} µm  ->  BBA = {row['BBA']:.4f}  "
            f"(spectral range [{alb.min():.4f}, {alb.max():.4f}])"
        )

    if args.plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        for _, row in df_spec.iterrows():
            ax.plot(wavelengths, row["albedo"], label=f"r = {int(row['rds'])} µm")
        ax.set_xlabel("Wavelength (µm)")
        ax.set_ylabel("Albedo")
        ax.set_xlim(0.2, 2.5)
        ax.set_ylim(0, 1.05)
        ax.set_title("Spectral albedo for different grain radii")
        ax.legend()
        fig.tight_layout()

    # Show all figures at once
    if args.plot:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()
