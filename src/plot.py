from classes import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def setup_axes(plot_config):

    rc = {
        "figure.figsize": (8, 6),
        "axes.facecolor": plot_config.facecolor,
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
    plt.style.use("seaborn")
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
        plt.savefig(str(savepath + "spectral_albedo.png"))

    if plot_config.show:
        plt.show()

    return
