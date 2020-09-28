import numpy as _np
from numpy import pi


def spheres_and_cylinders(target,
                          pore_diameter='pore.diameter',
                          throat_diameter='throat.diameter',
                          throat_length=None,
                          conduit_lengths=None,
                          return_elements=False):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders

    Parameter
    ---------
    target: OpenPNM object

    Notes
    -----
    The diffusive shape coefficient is the geometrical part of the pre-factor
    in Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    # Getting areas
    A1 = (pi/4*network[pore_diameter]**2)[cn[:, 0]]
    A2 = (pi/4*network[pore_diameter]**2)[cn[:, 1]]
    At = (pi/4*network[throat_diameter]**2)[throats]
    if conduit_lengths is not None:
        vals = target[conduit_lengths]
        L1, L2, Lt = vals['pore1'], vals['pore2'], vals['throat']
    else:
        a = target[throat_diameter][throats]
        r = target[pore_diameter][cn]
        theta = _np.arcsin(_np.atleast_2d(a).T/r)
        L1, L2 = (r*_np.cos(theta)).T
        if throat_length is not None:
            Lt = target[throat_length][throats]
        else:
            C1 = network['pore.coords'][cn[:, 0]]
            C2 = network['pore.coords'][cn[:, 1]]
            L = _np.sqrt(_np.sum((C1 - C2)**2, axis=1))
            Lt = L - L1 - L2
    # Find g for half of pore 1, the throat, and half of pore 2
    g1, g2, gt = A1/L1, A2/L2, At/Lt
    # Apply shape factors and calculate the final conductance
    if return_elements:
        vals = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        vals = (1/gt + 1/g1 + 1/g2)**(-1)
    return vals


def cylinders_in_series(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        n_cylinders=5,
                        throat_length=None,
                        return_elements=False):
    r"""

    """
    network = target.network


def pyramids_and_cuboids(target,
                         pore_diameter='pore.diameter',
                         throat_diameter='throat.diameter',
                         throat_length=None,
                         return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2


def cones_and_cylinders(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        throat_length='throat.length',
                        return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2
    Lt = target[throat_length]
    L1 = Lt - R1
    L2 = Lt - R2

    alpha1 = (R1 - Rt)/L1
    beta1 = 1 / (1/(Rt**3) - 1/(R1**3))
    alpha2 = (R2-Rt)/L2
    beta2 = 1 / (1/(Rt**3) - 1/(R2**3))
    g1 = (3*alpha1*pi/8) * beta1
    g2 = (3*alpha2*pi/8) * beta2
    gt = pi*Rt**4/(8*Lt)

    if return_elements:
        g = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        g = (1/g1 + 1/gt + 1/g2)**-1
    return g
