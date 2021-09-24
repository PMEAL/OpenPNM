import numpy as np
from openpnm.topotools import generate_base_points


def trim(network, pores=None, throats=None):
    r"""
    Remove given pores or throats from a network

    Parameters
    ----------
    network : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns', describing
        the network
    pores : array_like
        The list of pores to trim.
    throats : array_like
        The list of throats to trim

    Returns
    -------
    network : dict
        The ``network`` object with the specified pores/throats removed

    Notes
    -----
    Only one of ``pores`` or ``throats`` can be given.  To trim both types,
    call the fuction twice. This function renumbers the ``'throat.conns'``
    array when the pores being pointed to are removed.

    """
    if (pores is None) and (throats is None):
        raise Exception('Cannot trim pores and throats at the same time')
    if throats is not None:
        throats = np.atleast_1d(throats)
        keep = np.ones(network['throat.conns'].shape[0], dtype=bool)
        keep[throats] = False
        for item in network.keys():
            if item.startswith('throat'):
                network[item] = network[item][keep]
    elif pores is not None:
        pores = np.atleast_1d(pores)
        if pores.dtype == bool:
            pores = np.where(pores)[0]
        keep = np.ones(network['pore.coords'].shape[0], dtype=bool)
        keep[pores] = False
        for item in network.keys():
            if item.startswith('pore'):
                network[item] = network[item][keep]
        # Remove throats
        throats = np.any(np.isin(network['throat.conns'], pores), axis=1)
        network = trim(network, throats=throats)
        # Renumber throat conns
        remapping = np.cumsum(keep) - 1
        network['throat.conns'] = remapping[network['throat.conns']]
    return network


def join(net1, net2, L_max=0.99):
    r"""
    Joins two networks together topologically including new connections

    Parameters
    ----------
    net1 : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns'.
    net2 : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns'
    L_max : float
        The maximum distance between pores below which they are called
        neighbors

    Returns
    -------
    network : dict
        A dictionary containing 'pore.coords' with pores from both ``net1`` and
        ``net2``, and ``throat.conns`` with original connections plus new ones
        found during the join process.

    Notes
    -----
    This function uses ``scipy.spatial.KDTree``.

    """
    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1['pore.coords'])
    t2 = KDTree(net2['pore.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1['pore.coords'].shape[0]
    Np2 = net2['pore.coords'].shape[0]
    net3['pore.coords'] = np.vstack((net1.pop('pore.coords'),
                                     net2.pop('pore.coords')))
    net3['throat.conns'] = np.vstack((net1.pop('throat.conns'),
                                      net2.pop('throat.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3['throat.conns'] = np.vstack((net3.pop('throat.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, None)
        temp2 = net2.pop(item, None)
        if temp1 is None:
            temp1 = np.zeros(Np1, dtype=temp2.dtype)
        else:
            temp2 = np.zeros(Np2, dtype=temp1.dtype)
        net3[item] = np.hstack((temp1, temp2)).T
    return net3


def parse_points(shape, points):
    # Deal with input arguments
    if isinstance(points, int):
        points = generate_base_points(num_points=points,
                                      domain_size=shape,
                                      reflect=True)
    else:
        # Should we check to ensure that points are reflected?
        points = np.array(points)
    # Deal with points that are only 2D...they break Delaunay and Voronoi
    if points.shape[1] == 3 and len(np.unique(points[:, 2])) == 1:
        points = points[:, :2]

    return points


def get_spacing(network):
    r"""
    Determine spacing of a cubic network

    Parameters
    ----------
    network : dictionary
        A network dictionary containing 'pore.coords' and 'throat.conns'

    Returns
    -------
    spacing : ndarray
        An array containing the spacing between pores in each direction.

    Notes
    -----
    This function only works on simple cubic networks with no boundary pores.
    If a unique spacing cannot be found in each direction, and/or the throats
    are not all oriented perpendicularly, exceptions will be raised.

    """
    from openpnm.topotools import dimensionality
    # Find Network spacing
    P12 = network["throat.conns"]
    C12 = network["pore.coords"][P12]
    mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
    unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
    spacing = [0, 0, 0]
    dims = dimensionality(network)
    # Ensure vectors point in n-dims unique directions
    c = {tuple(row): 1 for row in unit_vec}
    mag = np.atleast_1d(mag.squeeze()).astype(float)
    if len(c.keys()) > sum(dims):
        raise Exception(
            "Spacing is undefined when throats point in more directions"
            " than network has dimensions."
        )
    for ax in [0, 1, 2]:
        if dims[ax]:
            inds = np.where(unit_vec[:, ax] == unit_vec[:, ax].max())[0]
            temp = np.unique(mag[inds])
            if not np.allclose(temp, temp[0]):
                raise Exception("A unique value of spacing could not be found.")
            spacing[ax] = temp[0]
    return np.array(spacing)


def get_shape(network):
    L = np.ptp(network["pore.coords"], axis=0)
    mask = L.astype(bool)
    S = get_spacing(network)
    shape = np.array([1, 1, 1], int)
    shape[mask] = L[mask] / S[mask] + 1
    return shape


def add_all_label(network):
    network['pore.all'] = np.ones(network['pore.coords'].shape[0], dtype=bool)
    network['throat.all'] = np.ones(network['throat.conns'].shape[0], dtype=bool)
    return network


def label_surface_pores(network):
    r"""
    """
    hits = np.zeros_like(network.Ps, dtype=bool)
    dims = topotools.dimensionality(network)
    mn = np.amin(network["pore.coords"], axis=0)
    mx = np.amax(network["pore.coords"], axis=0)
    for ax in [0, 1, 2]:
        if dims[ax]:
            hits += network["pore.coords"][:, ax] <= mn[ax]
            hits += network["pore.coords"][:, ax] >= mx[ax]
    network["pore.surface"] = hits
    return network
