from numpy import pi as _pi


def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    Returns
    -------
    value : NumPy ndarray
        Array containing pore area values.

    """
    diams = target[pore_diameter]
    value = _pi/4*(diams)**2
    return value


def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    Returns
    -------
    value : NumPy ndarray
        Array containing pore area values.

    """
    diams = target[pore_diameter]
    value = diams**2
    return value


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a circle

    Only suitable for true 2D simulations

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    """
    return target[pore_diameter]


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a square

    Only suitable for true 2D simulations

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.


    """
    return target[pore_diameter]
