import scipy as _sp


def unit_vector(target):
    r"""
    Calculates unit vector between each pair of pore centers.  In most
    networks this will correspond to the orientation of the throat connecting
    the pore-pairs.
    """
    network = target.project.network
    P12 = network['throat.conns']
    vec = _sp.diff(network['pore.coords'][P12], axis=1).squeeze()
    mag = _sp.atleast_2d(_sp.sqrt(_sp.sum(_sp.square(vec), axis=1))).T
    vals = vec/mag
    return vals[network.throats(target.name)]


def centroid(target, throat_endpoints='throat.endpoint'):
    r"""
    Find the centroid of each throat given the endpoint coordinates.

    Parameters
    ----------
    throat_endpoints : string
        An Nt x 2 x 3 array containing [[X1, Y1, Z1], [X2, Y2, Z2]] pairs on
        each row indicating the start and end points of a throat.  These
        should be arranged in the same order as 'throat.conns', so if pore 1
        is connected to pore 2, then this array should have the throat
        endpoint associated with pore 1 first, and with pore 2 second.
    """
    vals = _sp.mean(target['throat.endpoints'], axis=1)
    return vals
