r"""
===============================================================================
Submodule -- pore_surface_area
===============================================================================

"""
import scipy as _sp


def sphere(geometry, network,
           pore_diameter='pore.diameter',
           throat_area='throat.area',
           **kwargs):
    r"""
    Calculates internal surface area of pore bodies assuming they are spherical.
    then subtracts the area of the neighboring throats in a crude way, by
    simply considering the throat cross-sectional area, thus not accounting
    for the actual curvature of the intersection.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    network : OpenPNM Network Object
        The Network object associated with the Geometry.  This is needed to
        provide some topological information such as throat connections, and
        neighboring pores.

    pore_diameter : string
        The dictionary key to the pore diameter array.

    throat_area : string
        The dictioanry key to the throat area array.  Throat areas are needed
        since their insection with the pore are removed from the computation.

    """
    R = geometry[pore_diameter]/2
    Asurf = 4*_sp.constants.pi*R**2
    Tn = network.find_neighbor_throats(pores=geometry.Ps, flatten=False)
    Tsurf = _sp.array([_sp.sum(network[throat_area][Ts]) for Ts in Tn])
    value = Asurf - Tsurf
    return value


def cube(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate internal surface area for a cubic pore
    """
    D = geometry[pore_diameter]
    value = 6*D**2
    return value
