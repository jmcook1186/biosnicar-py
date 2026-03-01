"""Shared smoothing utility for radiative transfer solvers."""

from scipy.signal import savgol_filter


def apply_smoothing_function(albedo, model_config):
    """Applies Savitsky-Golay smoothing function to albedo, if toggled.

    Args:
        albedo: array of albedo values, likely passed as outputs.albedo
        model_config: instance of ModelConfig

    Returns:
        albedo: smoothed array of albedo values
    """
    return savgol_filter(albedo, model_config.window_size, model_config.poly_order)
