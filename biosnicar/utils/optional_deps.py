def import_miepython():
    """Import and return the miepython module. Intall hint on fail."""
    try:
        import miepython
    except ImportError as err:
        raise ImportError(
            "miepython is needed to generate SSOPs but is not installed. Install with `pip install biosnicar[mie]`"
            "or pip install miepython. It is not needed for the model itself, which reads precomputed files." 
        ) from err
    return miepython
