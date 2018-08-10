from numpy import pi as _pi
import numpy as _np


def spherical_pores(target, throat_area='throat.area',
                    pore_diameter='pore.diameter',
                    throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate equivalent cross-sectional area for a conduit consisting of two
    spherical pores and a constant cross-section throat. This area can be later
    used to calculate hydraulic or diffusive conductance of the conduit.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    d1 = network[pore_diameter][cn[:, 0]]
    d2 = network[pore_diameter][cn[:, 1]]
    L1 = target[throat_conduit_lengths+'.pore1']
    L2 = target[throat_conduit_lengths+'.pore2']
    return {'pore1': d1*L1*_pi / (2*_np.arctanh(2*L1/d1)),
            'throat': target[throat_area],
            'pore2': d2*L2*_pi / (2*_np.arctanh(2*L2/d2))}


def truncated_pyramid(target, throat_area='throat.area',
                      pore_diameter='pore.diameter'):
    r"""
    Calculate equivalent cross-sectional area for a conduit consisting of two
    truncated pyramid pores and a constant cross-section throat. This area can
    be later used to calculate hydraulic or diffusive conductance of the
    conduit.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    d1 = network[pore_diameter][cn[:, 0]]
    d2 = network[pore_diameter][cn[:, 1]]
    td = _np.sqrt(target[throat_area])
    return {'pore1': d1*td,
            'throat': target[throat_area],
            'pore2': d2*td}


def circular_pores(target, throat_area='throat.area',
                   pore_diameter='pore.diameter',
                   throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate equivalent cross-sectional area for a conduit consisting of two
    spherical pores and a constant cross-section throat. This area can be later
    used to calculate hydraulic or diffusive conductance of the conduit.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    d1 = network[pore_diameter][cn[:, 0]]
    d2 = network[pore_diameter][cn[:, 1]]
    L1 = target[throat_conduit_lengths+'.pore1']
    L2 = target[throat_conduit_lengths+'.pore2']
    return {'pore1': 2*L1 / (_np.arctan(2*L1/_np.sqrt(d1**2 - 4*L1**2))),
            'throat': target[throat_area],
            'pore2': 2*L2 / (_np.arctan(2*L2/_np.sqrt(d2**2 - 4*L2**2)))}


def boundary(target, pore_area='pore.area', throat_area='throat.area',
             label='pore.internal'):
    r"""
    A method to be applied to boundary throats that copies the values from the
    internal pore, if throat area values already exist use those
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    conns = network['throat.conns'][throats]
    labelled_ps = network[label][conns]
    p1 = conns[:, 0]
    p2 = conns[:, 1]
    # copy p index if label is false
    # label must be true for at least one pore
    p1[~labelled_ps[:, 0]] = p2[~labelled_ps[:, 0]]
    p2[~labelled_ps[:, 1]] = p1[~labelled_ps[:, 1]]
    ea_p1 = network[pore_area][p1]
    ea_p2 = network[pore_area][p2]
    if throat_area in target.keys():
        ea_ts = target[throat_area]
    else:
        ea_ts = ea_p1
    return {'pore1': ea_p1, 'throat': ea_ts, 'pore2': ea_p2}
