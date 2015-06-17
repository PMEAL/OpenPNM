r"""
===============================================================================
pore_topology -- functions for monitoring and adjusting topology
===============================================================================

"""
import scipy as _sp


def get_subscripts(network, shape, **kwargs):
    r"""
    Return the 3D subscripts (i,j,k) into the cubic network

    Parameters
    ----------
    shape : list
        The (i,j,k) shape of the network in number of pores in each direction

    """
    if network.num_pores('internal') != _sp.prod(shape):
        print('Supplied shape does not match Network size, cannot proceed')
    else:
        template = _sp.atleast_3d(_sp.empty(shape))
        a = _sp.indices(_sp.shape(template))
        i = a[0].flatten()
        j = a[1].flatten()
        k = a[2].flatten()
        ind = _sp.vstack((i, j, k)).T
        vals = _sp.ones((network.Np, 3))*_sp.nan
        vals[network.pores('internal')] = ind
        return vals


def adjust_spacing(network, new_spacing, **kwargs):
    r"""
    Adjust the the pore-to-pore lattice spacing on a cubic network

    Parameters
    ----------
    new_spacing : float
        The new lattice spacing to apply

    Notes
    -----
    At present this method only applies a uniform spacing in all directions.
    This is a limiation of OpenPNM Cubic Networks in general, and not of the
    method.
    """
    coords = network['pore.coords']
    try:
        spacing = network._spacing
        coords = coords/spacing*new_spacing
        network._spacing = new_spacing
    except:
        pass
    return coords


def reduce_coordination(network, z, mode='random', **kwargs):
    r"""
    Reduce the coordination number to the specified z value

    Parameters
    ----------
    z : int
        The coordination number or number of throats connected a pore

    mode : string, optional
        Controls the logic used to trim connections.  Options are:

        - 'random': (default) Throats will be randomly removed to achieve a
                    coordination of z
        - 'max': All pores will be adjusted to have a maximum coordination of z
                 (not implemented yet)

    Returns
    -------
    A label array indicating which throats should be trimmed to achieve desired
    coordination.

    Notes
    -----
    Pores with only 1 throat will be ignored in all calculations since these
    are generally boundary pores.

    """
    T_trim = ~network['throat.all']
    T_nums = network.num_neighbors(network.pores())
    # Find protected throats
    T_keep = network.find_neighbor_throats(pores=(T_nums == 1))
    if mode == 'random':
        z_ave = _sp.average(T_nums[T_nums > 1])
        f_trim = (z_ave - z)/z_ave
        T_trim = _sp.rand(network.Nt) < f_trim
        T_trim = T_trim*(~network.tomask(throats=T_keep))
    if mode == 'max':
        pass
    return T_trim
