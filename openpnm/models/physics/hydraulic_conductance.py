r"""

.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille

"""

from .misc import generic_conductance
import numpy as np


def hagen_poiseuille(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_shape_factors : string
        Dictionary key of the conduit FLOW shape factor values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

    """
    return generic_conductance(target=target, transport_type='flow',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)


def valvatne_blunt(target,
                   pore_viscosity='pore.viscosity',
                   throat_viscosity='throat.viscosity',
                   pore_shape_factor='pore.shape_factor',
                   throat_shape_factor='throat.shape_factor',
                   pore_area='pore.area',
                   throat_area='throat.area',
                   conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the single phase hydraulic conductance of conduits in network,
    where a conduit is ( 1/2 pore - full throat - 1/2 pore ) according to [1].
    Function has been adapted for use with the Statoil imported networks and
    makes use of the shape factor in these networks to apply Hagen-Poiseuille
    flow for conduits of different shape classes: Triangular, Square and
    Circular [2].

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    pore_shape_factor : string
        Dictionary key of the pore geometric shape factor values

    throat_shape_factor : string
        Dictionary key of the throat geometric shape factor values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_lengths : string
        Dictionary key of the throat conduit lengths

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    References
    ----------
    [1] Valvatne, Per H., and Martin J. Blunt. "Predictive pore‐scale modeling
    of two‐phase flow in mixed wet media." Water Resources Research 40,
    no. 7 (2004).
    [2] Patzek, T. W., and D. B. Silin (2001), Shape factor and hydraulic
    conductance in noncircular capillaries I. One-phase creeping flow,
    J. Colloid Interface Sci., 236, 295–304.
    """
    network = target.project.network
    mu_p = target[pore_viscosity]
    try:
        mu_t = target[throat_viscosity]
    except KeyError:
        mu_t = target.interpolate_data(pore_viscosity)
    # Throat Portion
    Gt = network[throat_shape_factor]
    tri = Gt <= np.sqrt(3)/36.0
    circ = Gt >= 0.07
    square = ~(tri | circ)
    ntri = np.sum(tri)
    nsquare = np.sum(square)
    ncirc = np.sum(circ)
    kt = np.ones_like(Gt)
    kt[tri] = 3.0/5.0
    kt[square] = 0.5623
    kt[circ] = 0.5
    # Pore Portions
    Gp = network[pore_shape_factor]
    tri = Gp <= np.sqrt(3)/36.0
    circ = Gp >= 0.07
    square = ~(tri | circ)
    ntri += np.sum(tri)
    nsquare += np.sum(square)
    ncirc += np.sum(circ)
    kp = np.ones_like(Gp)
    kp[tri] = 3.0/5.0
    kp[square] = 0.5623
    kp[circ] = 0.5
    gp = kp*(network[pore_area]**2)*Gp/mu_p
    gt = kt*(network[throat_area]**2)*Gt/mu_t
    conns = network['throat.conns']
    l1 = network[conduit_lengths + '.pore1']
    lt = network[conduit_lengths + '.throat']
    l2 = network[conduit_lengths + '.pore2']
    # Resistors in Series
    value = (l1/gp[conns[:, 0]] +
             lt/gt +
             l2/gp[conns[:, 1]])
    return 1/value
