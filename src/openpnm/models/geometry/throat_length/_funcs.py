from openpnm.models.geometry import _geodocs


__all__ = ["spheres_and_cylinders",
           "circles_and_rectangles",
           "cones_and_cylinders",
           "intersecting_cones",
           "hybrid_cones_and_cylinders",
           "trapezoids_and_rectangles",
           "hybrid_trapezoids_and_rectangles",
           "intersecting_trapezoids",
           "pyramids_and_cuboids",
           "intersecting_pyramids",
           "hybrid_pyramids_and_cuboids",
           "cubes_and_cuboids",
           "squares_and_rectangles"]


@_geodocs
def spheres_and_cylinders(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are spheres and throats are
    cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------
    lengths : ndarray
        A numpy ndarray containing throat length values

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.spheres_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def circles_and_rectangles(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter',
):
    r"""
    Finds throat length assuming pores are circles and throats are
    rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.circles_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def cones_and_cylinders(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cones and throats are
    cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def intersecting_cones(
    network,
    pore_coords="pore.coords",
    throat_coords="throat.coords"
):
    r"""
    Finds throat length assuming pores are intersecting cones.

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.intersecting_cones(
        network=network,
        pore_coords=pore_coords,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def hybrid_cones_and_cylinders(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Finds throat length assuming pores are cones and throats are
    cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.hybrid_cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def trapezoids_and_rectangles(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are trapezoids and throats are
    rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.trapezoids_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def intersecting_trapezoids(
    network,
    pore_coords="pore.coords",
    throat_coords="throat.coords"
):
    r"""
    Finds throat length assuming pores are intersecting trapezoids.

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.intersecting_trapezoids(
        network=network,
        pore_coords=pore_coords,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def hybrid_trapezoids_and_rectangles(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Finds throat length assuming pores are trapezoids and throats are
    rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.hybrid_trapezoids_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def pyramids_and_cuboids(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are pyramids and throats are
    cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.pyramids_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def intersecting_pyramids(
    network,
    pore_coords="pore.coords",
    throat_coords="throat.coords"
):
    r"""
    Finds throat length assuming pores are intersecting pyramids.

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.intersecting_pyramids(
        network=network,
        pore_coords=pore_coords,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def hybrid_pyramids_and_cuboids(
    network,
    pore_diameter='pore.diameter',
    throat_coords="throat.coords",
):
    r"""
    Finds throat length assuming pores are pyramids and throats are
    cuboids.

    Parameters
    ----------
    %(network)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.hybrid_pyramids_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    )
    return out[:, 1]


@_geodocs
def cubes_and_cuboids(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cubes and throats are
    cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cubes_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]


@_geodocs
def squares_and_rectangles(
    network,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are squares and throats are
    rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.squares_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )
    return out[:, 1]
