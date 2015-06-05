r"""
===============================================================================
pore_area -- Models for cross-sectional area of a pore body
===============================================================================

"""
import scipy as _sp


def spherical(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[10,10,10])
    >>> geo = OpenPNM.Geometry.GenericGeometry(network=pn,
                                               pores=pn.Ps,
                                               throats=pn.Ts)
    >>> geo['pore.diameter'] = sp.rand(geo.Np)
    >>> list(geo.props())  # Check that the seed_values are present
    ['pore.diameter']
    >>> geo.models.add(propname = 'pore.area',
                       model = OpenPNM.Geometry.models.pore_area.sphere,
                       pore_diameter = 'pore.diameter')
    >>> sorted(list(geo.models))  # Check that the model is present
    ['pore.area']
    >>> sorted(list(geo.props()))  # Check that the numerical values are there
    ['pore.diameter', 'pore.area']

    """
    diams = geometry[pore_diameter]
    value = _sp.constants.pi/4*(diams)**2
    return value


def cubic(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[10,10,10])
    >>> geo = OpenPNM.Geometry.GenericGeometry(network=pn,
                                               pores=pn.Ps,
                                               throats=pn.Ts)
    >>> geo['pore.diameter'] = sp.rand(geo.Np)
    >>> list(geo.props())  # Check that the seed_values are present
    ['pore.diameter']
    >>> geo.models.add(propname = 'pore.area',
                       model = OpenPNM.Geometry.models.pore_area.cube,
                       pore_diameter = 'pore.diameter')
    >>> sorted(list(geo.models))  # Check that the model is present
    ['pore.area']
    >>> sorted(list(geo.props()))  # Check that the numerical values are there
    ['pore.diameter', 'pore.area']

    """
    diams = geometry[pore_diameter]
    value = diams**2
    return value
