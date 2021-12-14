"""
Tools
-----

"""
import numpy as np
from openpnm.topotools import generate_base_points, dimensionality, isoutside


def trim(network, vert_ids=None, edge_ids=None):
    r"""
    Remove given pores or throats from a network

    Parameters
    ----------
    network : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns', describing
        the network
    vert_ids : array_like
        The list of vertcies to trim.
    edge_ids : array_like
        The list of edges to trim

    Returns
    -------
    network : dict
        The ``network`` object with the specified vertices/edges removed

    Notes
    -----
    Only one of ``vert_ids`` or ``edge_ids`` can be given.  To trim both types,
    call the fuction twice. This function renumbers the ``'edge.conns'``
    array when the vertices being pointed to are removed.

    """
    if (vert_ids is not None) and (edge_ids is not None):
        raise Exception('Cannot trim pores and throats at the same time')
    if edge_ids is not None:
        edge_ids = np.atleast_1d(edge_ids)
        keep = np.ones(network['edge.conns'].shape[0], dtype=bool)
        keep[edge_ids] = False
        for item in network.keys():
            if item.startswith('edge'):
                network[item] = network[item][keep]
    elif vert_ids is not None:
        vert_ids = np.atleast_1d(vert_ids)
        if vert_ids.dtype == bool:
            vert_ids = np.where(vert_ids)[0]
        keep = np.ones(network['vert.coords'].shape[0], dtype=bool)
        keep[vert_ids] = False
        for item in network.keys():
            if item.startswith('vert'):
                network[item] = network[item][keep]
        # Remove edges
        edges = np.any(np.isin(network['edge.conns'], vert_ids), axis=1)
        network = trim(network, edge_ids=edges)
        # Renumber throat conns
        remapping = np.cumsum(keep) - 1
        network['edge.conns'] = remapping[network['edge.conns']]
    return network


def join(net1, net2, L_max=0.99):
    r"""
    Joins two networks together topologically including new connections

    Parameters
    ----------
    net1 : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns'.
    net2 : dictionary
        A dictionary containing 'vert.coords' and 'edge.conns'
    L_max : float
        The maximum distance between vertices below which they are called
        neighbors

    Returns
    -------
    network : dict
        A dictionary containing 'vert.coords' vertices from both ``net1`` and
        ``net2``, and ``edge.conns`` with original connections plus new ones
        found during the join process.

    Notes
    -----
    This function uses ``scipy.spatial.KDTree``.

    """
    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1['vert.coords'])
    t2 = KDTree(net2['vert.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1['vert.coords'].shape[0]
    Np2 = net2['vert.coords'].shape[0]
    net3['vert.coords'] = np.vstack((net1.pop('vert.coords'),
                                     net2.pop('vert.coords')))
    net3['edge.conns'] = np.vstack((net1.pop('edge.conns'),
                                    net2.pop('edge.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3['edge.conns'] = np.vstack((net3.pop('edge.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, np.zeros(Np1)*np.nan)
        temp2 = net2.pop(item, np.zeros(Np2)*np.nan)
        net3[item] = np.concatenate((temp1, temp2), axis=0)
    return net3


def parse_points(shape, points, reflect=False):
    # Deal with input arguments
    if isinstance(points, int):
        points = generate_base_points(num_points=points,
                                      domain_size=shape,
                                      reflect=reflect)
    else:
        # Should we check to ensure that points are reflected?
        points = np.array(points)
    return points


def get_spacing(network):
    r"""
    Determine spacing of a cubic network

    Parameters
    ----------
    network : dictionary
        A network dictionary containing 'vert.coords' and 'edge.conns'

    Returns
    -------
    spacing : ndarray
        An array containing the spacing between vertices in each direction

    Notes
    -----
    This function only works on simple cubic networks with no boundary
    vertices. If a unique spacing cannot be found in each direction,
    and/or the edges are not all oriented perpendicularly, exceptions
    will be raised.

    """
    from openpnm.topotools import dimensionality
    # Find Network spacing
    P12 = network["edge.conns"]
    C12 = network["vert.coords"][P12]
    mag = np.linalg.norm(np.diff(C12, axis=1), axis=2)
    unit_vec = np.around(np.squeeze(np.diff(C12, axis=1)) / mag, decimals=14)
    spacing = [0, 0, 0]
    dims = dimensionality(coords=network['vert.coords'])
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
    L = np.ptp(network["vert.coords"], axis=0)
    mask = L.astype(bool)
    S = get_spacing(network)
    shape = np.array([1, 1, 1], int)
    shape[mask] = L[mask] / S[mask] + 1
    return shape


def add_all_label(network):
    network['pore.all'] = np.ones(network['pore.coords'].shape[0], dtype=bool)
    network['throat.all'] = np.ones(network['throat.conns'].shape[0], dtype=bool)
    return network


def label_surface_nodes(network):
    r"""
    """
    hits = np.zeros_like(network.Ps, dtype=bool)
    dims = dimensionality(network)
    mn = np.amin(network["vert.coords"], axis=0)
    mx = np.amax(network["vert.coords"], axis=0)
    for ax in np.where(dims)[0]:
        if dims[ax]:
            hits += network["vert.coords"][:, ax] <= mn[ax]
            hits += network["vert.coords"][:, ax] >= mx[ax]
    network["vert.surface"] = hits
    return network


def label_faces(network, threshold=0.05):
    r"""
    Label the vertices sitting on the faces of the domain in accordance with
    the conventions used for cubic etc.
    """
    dims = dimensionality(network)
    coords = np.around(network['vert.coords'], decimals=10)
    min_labels = ['front', 'left', 'bottom']
    max_labels = ['back', 'right', 'top']
    min_coords = np.amin(coords, axis=0)
    max_coords = np.amax(coords, axis=0)
    for ax in np.where(dims)[0]:
        network['vert.' + min_labels[ax]] = coords[:, ax] <= threshold*min_coords[ax]
        network['vert.' + max_labels[ax]] = coords[:, ax] >= (1-threshold)*max_coords[ax]
    return network


def crop(network, shape, mode='full'):
    r"""
    Removes vertices that lie outside the specifed shape

    Parameters
    ----------
    network : dict
        Dictionary containing ``vert.coords`` and ``edge.conns`` arrays
    shape : array_like
        The [x, y, z] shape of the domain beyond which trimming should be
        applied
    mode : str
        Controls how vertices to be trimmed is determined. Options are:

            * 'full':
                Any vertices lying outside the domain are trimmed
            * 'mixed'
                Vertices with at least one neighbor lying inside the domain
                are kept
    """
    Pdrop = isoutside(network['vert.coords'], shape=shape, thresh=0)
    if mode == 'full':
        network = trim(network=network, pores=np.where(Pdrop)[0])
    elif mode == 'mixed':
        # Find throats connecting internal to external pores
        Ts = np.sum(Pdrop[network['edge.conns']], axis=1) == 1
        # Keep the pores on the ends of these throats
        Pkeep = np.unique(network['edge.conns'][Ts])
        Ps = np.array(list(set(np.where(Pdrop)[0]).difference(set(Pkeep)))).astype(int)
        network = trim(network=network, pores=Ps)
        # Remove throats between these surviving external pores
        Pdrop = isoutside(network['vert.coords'], shape=shape)
        Ts = np.all(Pdrop[network['edge.conns']], axis=1)
        network = trim(network=network, throats=Ts)
        # Lastly label the surviving pores as outside for further processing
        network['vert.outside'] = np.zeros(network['vert.coords'].shape[0], dtype=bool)
        network['vert.outside'][Pdrop] = True
    return network
