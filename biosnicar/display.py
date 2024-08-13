#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def setup_axes(plot_config):
    rc = {
        "figure.figsize": (8, 6),
        "axes.facecolor": str(plot_config.facecolor),
        "axes.grid": plot_config.grid,
        "grid.color": str(plot_config.grid_color),
        "xtick.major.width": plot_config.xtick_width,
        "xtick.major.size": plot_config.xtick_size,
        "ytick.major.width": plot_config.ytick_width,
        "ytick.major.size": plot_config.ytick_size,
        "axes.linewidth": plot_config.linewidth,
        "font.size": plot_config.fontsize,
        "xtick.bottom": plot_config.xtick_btm,
        "ytick.left": plot_config.ytick_left,
    }

    return rc


def plot_albedo(plot_config, model_config, albedo):
    rc = setup_axes(plot_config)
    sns.set_style("white")
    plt.rcParams.update(rc)
    wvl = model_config.wavelengths

    plt.plot(wvl, albedo), plt.ylabel("Albedo", fontsize=18), plt.xlabel(
        "Wavelength (µm)", fontsize=18
    ), plt.xlim(0.3, 2.5), plt.ylim(0, 1), plt.xticks(fontsize=15), plt.yticks(
        fontsize=15
    ), plt.axvline(
        x=0.68, color="g", linestyle="dashed"
    )

    if plot_config.save:
        plt.savefig(str(model_config.savefigpath + "spectral_albedo.png"))

    plt.show()

    return


def display_out_data(outputs):
    print("\n** OUTPUT DATA **")
    print("Broadband albedo: ", np.round(outputs.BBA, 4))

    I2DBA, I3DBA, NDCI, MCI, II = calculate_band_ratios(outputs.albedo)

    print("\nBAND RATIO INDEX VALUES")
    print("2DBA Index: ", I2DBA)
    print("3DBA index: ", I3DBA)
    print("NDCI index: ", NDCI)
    print("MCI index: ", MCI)
    print("Impurity Index: ", II)

    return


def calculate_band_ratios(albedo):
    I2DBA = albedo[51] / albedo[46]
    I3DBA = (albedo[46] - albedo[50]) / albedo[55]
    NDCI = ((albedo[50] - albedo[48]) - (albedo[55] - albedo[48])) * (
        (albedo[50] - albedo[48]) / (albedo[55] - albedo[48])
    )
    MCI = (albedo[50] - albedo[46]) / (albedo[50] + albedo[46])
    II = np.log(albedo[36]) / np.log(albedo[66])

    return I2DBA, I3DBA, NDCI, MCI, II


if __name__ == "__main__":
    pass
