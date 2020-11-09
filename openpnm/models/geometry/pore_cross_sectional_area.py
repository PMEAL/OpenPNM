from numpy import pi as _pi


def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere.

    Parameters
    ----------
    target : GenericGeometry
        The Geometry object which this model is associated with. This
        controls the length of the calculated array, and also provides
        access to other necessary geometric properties.
    pore_diameter : str
        The dictionary key of the array on the Geometry object containing
        the pore diameter values necessary to find the area.

    Returns
    -------
    ndarray
        Array containing pore area values.

    """
    D = target[pore_diameter]
    return _pi/4 * D**2


def cone(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cone.

    Parameters
    ----------
    target : GenericGeometry
        The Geometry object which this model is associated with. This
        controls the length of the calculated array, and also provides
        access to other necessary geometric properties.
    pore_diameter : str
        The dictionary key of the array on the Geometry object containing
        the pore diameter values necessary to find the area.

    Returns
    -------
    ndarray
        Array containing pore area values.

    """
    D = target[pore_diameter]
    return _pi / 4 * D**2


def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    target : GenericGeometry
        The Geometry object which this model is associated with. This
        controls the length of the calculated array, and also provides
        access to other necessary geometric properties.
    pore_diameter : str
        The dictionary key of the array on the Geometry object containing
        the pore diameter values necessary to find the area.

    Returns
    -------
    ndarray
        Array containing pore area values.

    """
    D = target[pore_diameter]
    return D**2


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a circle.

    Parameters
    ----------
    target : GenericGeometry
        The Geometry object which this model is associated with. This
        controls the length of the calculated array, and also provides
        access to other necessary geometric properties.
    pore_diameter : str
        The dictionary key of the array on the Geometry object containing
        the pore diameter values necessary to find the area.

    Returns
    -------
    ndarray
        Array containing pore area values.

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return target[pore_diameter]


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a square.

    Parameters
    ----------
    target : GenericGeometry
        The Geometry object which this model is associated with. This
        controls the length of the calculated array, and also provides
        access to other necessary geometric properties.
    pore_diameter : str
        The dictionary key of the array on the Geometry object containing
        the pore diameter values necessary to find the area.

    Returns
    -------
    ndarray
        Array containing pore area values.

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return target[pore_diameter]
