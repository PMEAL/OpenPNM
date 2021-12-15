from openpnm.utils import Docorator


__all__ = ["spheres_and_cylinders",
           "circles_and_rectangles",
           "cones_and_cylinders",
           "trapezoids_and_rectangles",
           "pyramids_and_cuboids",
           "cubes_and_cuboids",
           "squares_and_rectangles"]
docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_length',
                     sections=['Parameters', 'Returns'])
def spheres_and_cylinders(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are spheres and throats are
    cylinders.

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Name of the dictionary key on ``target`` where the array containing
        pore diameter values is stored
    throat_diameter : str
        Name of the dictionary key on ``target`` where the array containing
        throat diameter values is stored

    Returns
    -------
    lengths : ndarray
        A numpy ndarray containing throat length values

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.spheres_and_cylinders(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def circles_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are circles and throats are
    rectangles.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns)

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.circles_and_rectangles(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def cones_and_cylinders(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cones and throats are
    cylinders.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cones_and_cylinders(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def trapezoids_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are trapezoids and throats are
    rectangles.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.trapezoids_and_rectangles(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def pyramids_and_cuboids(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are pyramids and throats are
    cuboids.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.pyramids_and_cuboids(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def cubes_and_cuboids(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cubes and throats are
    cuboids.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cubes_and_cuboids(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


def squares_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are squares and throats are
    rectangles.

    Parameters
    ----------
    %(models.geometry.throat_length.parameters)

    Returns
    -------
    %(models.geometry.throat_length.returns

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.squares_and_rectangles(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]
