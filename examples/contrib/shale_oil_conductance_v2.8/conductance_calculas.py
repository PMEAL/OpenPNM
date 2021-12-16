# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 10:13:57 2021

@author: xu kai
"""

import numpy as np

def slip_shale_conductance(target,
        negative_slippage = False,
        pore_area = "pore.cross_sectional_area",
        throat_area = "throat.cross_sectional_area",
        pore_area_a = "pore.cross_sectional_area_a",
        throat_area_a = "throat.cross_sectional_area_a",
        pore_viscosity = "pore.viscosity",
        throat_viscosity = "throat.viscosity",
        pore_viscosity_a = "pore.viscosity_a",
        throat_viscosity_a = "throat.viscosity_a",
        pore_slip_length = "pore.l_sd",
        throat_slip_length = "throat.l_sd",
        conduit_lengths = "throat.conduit_lengths",
        pore_shapefactor = "pore.shape_factor",
        throat_shapefactor = "throat.shape_factor"):

    """
    Calculate conductance for throats for liquid in shale
    noncircular nanopores, considering slip effect.
    positive slippage or negative slippage

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    negative_slippage: Bool
        if True, should consider the situation where negative slippage exists.
        default is False.

    pore_area : string
        The dictionary key of the array containing the pore area values.

    throat_area : string
        The dictionary key of the array containing the throat area values.

    pore_area_a : string
        The dictionary key of the array containing the pore adsorption
        layer area values of pores.

    throat_area_a : string
        The dictionary key of the array containing the throat adsorption
        layer area values of pores.

    pore_viscosity : string
        The dictionary key of the array containing the pore viscosity values.

    throat_viscosity : string
        The dictionary key of the array containing the throat viscosity values.

    pore_viscosity_a : string
        The dictionary key of the array containing the adsorption layer
        viscosity values of pores.

    throat_viscosity_a : string
        The dictionary key of the array containing the adsorption layer
        viscosity value of throats.

    pore_slip_length : string
        The dictionary key of the array containing the slip length in pores.

    throat_slip_length : string
        The dictionary key of the array containing the slip length in throats.

    conduit_lengths : string
        The dictionary key of the array containing the conduit length values.

    pore_shapefactor : string
        The dictionary key of the array containing the shape factor
        values of pores.
    throat_shapefactor : string
        The dictionary key of the array containing the shape factor
        values of throats.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat hydraulic conductance.

    References
    -----------
    Afsharpoor A., Javadpour F. Liquid slip flow in a network of shale
    noncircular nanopores[J]. Fuel, 2016, 180: 580-590.

    Yang Y., Wang K., Zhang L., et al. Pore-scale simulation of shale oil 
    flow based on pore network model[J]. Fuel, 2019, 251: 683-692. 
    """

    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    A1_a = network[pore_area_a][cn[:, 0]]
    At_a = network[throat_area_a][throats]
    A2_a = network[pore_area_a][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths][:, 0]
    Lt = network[conduit_lengths][:, 1]
    L2 = network[conduit_lengths][:, 2]

    # Getting shape factors
    G1 = network[pore_shapefactor][cn[:, 0]]
    Gt = network[throat_shapefactor][throats]
    G2 = network[pore_shapefactor][cn[:, 1]]

    # Getting dimensionless slip length
    L1_sd = phase[pore_slip_length][cn[:, 0]]
    Lt_sd = phase[throat_slip_length][throats]
    L2_sd = phase[pore_slip_length][cn[:, 1]]

    # Preallocating g
    g1, g2, gt = np.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = np.inf

    # bulk viscosity
    Dt = phase[throat_viscosity][throats]
    D1, D2 = phase[pore_viscosity][cn].T

    # adsorption layer viscosity
    Dt_a = phase[throat_viscosity_a][throats]
    D1_a, D2_a = phase[pore_viscosity_a][cn].T

    # equivalent viscosity
    mu_1 = (D1 * (A1 - A1_a) + D1_a * A1_a) / A1
    mu_2 = (D2 * (A2 - A2_a) + D2_a * A2_a) / A2
    mu_t = (Dt * (At - At_a) + Dt_a * At_a) / At

    # shape factor conditions
    c1, c2, ct = [G >= 0.04 for G in [G1, G2, Gt]]

    if negative_slippage == False:
        # conduit length and shape factor condition
        # when G >= 0.04
        t1 = np.logical_and(m1, c1)
        t2 = np.logical_and(m2, c2)
        tt = np.logical_and(mt, ct)

        fac1 = (-0.16 + 0.12*L1_sd + 6.4*G1 -5.5E-3*L1_sd**2 - 50*G1**2 + 1.7*L1_sd*G1)[t1]
        fac2 = (-0.16 + 0.12*L2_sd + 6.4*G2 -5.5E-3*L2_sd**2 - 50*G2**2 + 1.7*L2_sd*G2)[t2]
        fact = (-0.16 + 0.12*Lt_sd + 6.4*Gt -5.5E-3*Lt_sd**2 - 50*Gt**2 + 1.7*Lt_sd*Gt)[tt]

        # Find g for half of pore 1, throat, and half of pore 2
        g1[t1] = A1[t1] ** 2 * fac1 / (mu_1 * L1)[t1]
        g2[t2] = A2[t2] ** 2 * fac2 / (mu_2 * L2)[t2]
        gt[tt] = At[tt] ** 2 * fact / (mu_t * Lt)[tt]

        # G < 0.04
        t1 = np.logical_and(m1, ~c1)
        t2 = np.logical_and(m2, ~c2)
        tt = np.logical_and(mt, ~ct)

        fac1 = (-1.2E-2 + 5.7E-2*L1_sd + 2*G1 -5.2E-3*L1_sd**2 - 38*G1**2 + 3.2*L1_sd*G1)[t1]
        fac2 = (-1.2E-2 + 5.7E-2*L2_sd + 2*G2 -5.2E-3*L2_sd**2 - 38*G2**2 + 3.2*L2_sd*G2)[t2]
        fact = (-1.2E-2 + 5.7E-2*Lt_sd + 2*Gt -5.2E-3*Lt_sd**2 - 38*Gt**2 + 3.2*Lt_sd*Gt)[tt]

        # Find g for half of pore 1, throat, and half of pore 2
        g1[t1] = A1[t1] ** 2 * fac1 / (mu_1 * L1)[t1]
        g2[t2] = A2[t2] ** 2 * fac2 / (mu_2 * L2)[t2]
        gt[tt] = At[tt] ** 2 * fact / (mu_t * Lt)[tt]

    else:
        # conduit length and shape factor condition
        # when G >= 0.04
        t1 = np.logical_and(m1, c1)
        t2 = np.logical_and(m2, c2)
        tt = np.logical_and(mt, ct)

        fac1 = (-0.16 + 6.4*G1 - 50*G1**2)[t1]
        fac2 = (-0.16 + 6.4*G2 - 50*G2**2)[t2]
        fact = (-0.16 + 6.4*Gt - 50*Gt**2)[tt]

        # Find g for half of pore 1, throat, and half of pore 2
        g1[t1] = (A1 - A1_a)[t1] ** 2 * fac1 / (D1 * L1)[t1]
        g2[t2] = (A2 - A2_a)[t2] ** 2 * fac2 / (D2 * L2)[t2]
        gt[tt] = (At - At_a)[tt] ** 2 * fact / (Dt * Lt)[tt]

        # G < 0.04
        t1 = np.logical_and(m1, ~c1)
        t2 = np.logical_and(m2, ~c2)
        tt = np.logical_and(mt, ~ct)

        fac1 = (-1.2E-2 + 2*G1 - 38*G1**2)[t1]
        fac2 = (-1.2E-2 + 2*G2 - 38*G2**2)[t2]
        fact = (-1.2E-2 + 2*Gt - 38*Gt**2)[tt]

        g1[t1] = (A1 - A1_a)[t1] ** 2 * fac1 / (D1 * L1)[t1]
        g2[t2] = (A2 - A2_a)[t2] ** 2 * fac2 / (D2 * L2)[t2]
        gt[tt] = (At - At_a)[tt] ** 2 * fact / (Dt * Lt)[tt]

    # Apply shape factors and calculate the final conductance
    return (1/gt + 1/g1 + 1/g2) ** (-1)
