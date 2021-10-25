r"""
"""


def water(target):
    r"""
    Calculates the relative permittivity of saline water.

    Notes
    -----
    The relative dielectric constant or relative permittivity (sometimes simply
    called permittivity) of a medium is the ratio between the absolute
    permittivity to the vacuum permittivity. The vacuum permittivity or
    electric constant has a value of 8.854187817e-12 F/m.
    At 25°C, the relative permittivity of water is 78.303.

    Warnings
    --------
    The "ChargeConservation" algorithm (used by "IonicTransport" and coupled
    with "NernstPlanck") is implemented, when assumption of
    "charge_conservation" is "poisson", with the idea that the medium has a
    uniform dielectric constant. If the dielectric constant is not uniform,
    the permittivity should be removed from the source term and included into
    the "ionic_conductance" model.

    """
    value = 78.303
    return value
