import scipy as sp
import scipy.ndimage as spim
import scipy.sparse as sprs
from scipy.sparse import csgraph
from openpnm.core import logging, Workspace
from openpnm.utils.misc import PrintableDict
ws = Workspace()
logger = logging.getLogger()


def find_neighbor_sites(sites, am, flatten=True, exclude_input=True,
                        logic='union'):
    r"""
    Given a symmetric adjacency matrix, finds all sites that are connected
    to the input sites.

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    flatten : boolean (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.

    exclude_input : boolean (default is ``True``)
        If ``flatten`` is ``True``, then this flag will remove the input sites
        from the result, otherwise it's ignored.

    logic : string
        Specifies logic to apply to a flattened list.  Options are:

        **'union'** : (default) All neighbors of the input sites

        **'intersection'** : Only neighbors shared by all input sites

        **'exclusive_or'** : Only neighbors not shared by any input sites

    """
    if am.format != 'lil':
        am = am.tolil(copy=False)
    neighbors = [am.rows[i] for i in sp.array(sites, ndmin=1)]
    if flatten:
        neighbors.append(sites)
        neighbors = sp.hstack(neighbors)
        neighbors = apply_logic(neighbors=neighbors, logic=logic)
        if exclude_input:
            neighbors = neighbors[sp.in1d(neighbors, sites, invert=True,
                                          assume_unique=True)]
    return neighbors


def find_neighbor_bonds(sites, im, flatten=True, logic='union'):
    r"""
    Given an incidence matrix, finds all sites that are connected to the
    input sites.

    Parameters
    ----------
    im : scipy.sparse matrix
        The incidence matrix of the network.  Must be shaped as (N-sites,
        N-bonds), with non-zeros indicating which bonds are connected.

    flatten : boolean (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.

    logic : string
        Specifies logic to apply to a flattened list.  Options are:

        **'union'** : (default) All neighbors of the input sites

        **'intersection'** : Only neighbors shared by all input sites

        **'exclusive_or'** : Only neighbors not shared by any input sites

    """
    if im.shape[0] > im.shape[1]:
        print('Warning: Received matrix has more sites than bonds!')
    if im.format != 'lil':
        im = im.tolil(copy=False)
    neighbors = [im.rows[i] for i in sp.array(sites, ndmin=1)]
    if flatten:
        neighbors = sp.concatenate(neighbors)
        neighbors = apply_logic(neighbors=neighbors, logic=logic)
    return neighbors


def find_connected_sites(bonds, am, flatten=True, logic='union'):
    r"""
    Given an adjacency matrix, finds which sites are connected to the input
    bonds.

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    flatten : boolean (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.

    logic : string
        Specifies logic to apply to a flattened list.  Options are:

        **'union'** : (default) All neighbors of the input sites

        **'intersection'** : Only neighbors shared by all input sites

        **'exclusive_or'** : Only neighbors not shared by any input sites
    """
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    bonds = sp.array(bonds, ndmin=1)
    if flatten:
        neighbors = sp.hstack((am.row[bonds], am.col[bonds]))
        neighbors = apply_logic(neighbors=neighbors, logic=logic)
    else:
        neighbors = sp.vstack((am.row[bonds], am.col[bonds])).T
        neighbors = sp.array(neighbors, dtype=int)
    return neighbors


def find_connecting_bonds(sites, am):
    r"""
    Given pairs of sites, finds the bonds which connects each pair.

    Parameters
    ----------
    sites : array_like
        A 2-column vector containing pairs of site indices on each row.

    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    Returns
    -------
    Returns a list the same length as P1 (and P2) with the each element
    containing the throat number that connects the corresponding pores,
    or `None`` if pores are not connected.

    Notes
    -----
    The returned list can be converted to an ND-array, which will convert
    the ``None`` values to ``nan``.  These can then be found using
    ``scipy.isnan``.

    """
    if am.format != 'dok':
        am = am.todok(copy=False)
    sites = sp.array(sites, ndmin=2)
    z = tuple(zip(sites[:, 0], sites[:, 1]))
    neighbors = [am.get(z[i], None) for i in range(len(z))]
    return neighbors


def apply_logic(neighbors, logic):
    if neighbors.ndim > 1:
        neighbors = sp.hstack(neighbors)
    if logic == 'union':
        neighbors = sp.unique(neighbors)
    elif logic == 'exclusive_or':
        # neighbors that occur only once are NOT shared with other sites
        neighbors = sp.unique(sp.where(sp.bincount(neighbors) == 1)[0])
    elif logic == 'intersection':
        # neighbors that occur more than once ARE shared with other sites
        neighbors = sp.unique(sp.where(sp.bincount(neighbors) > 1)[0])
    else:
        raise Exception('Unsupported logic type: '+logic)
    return neighbors.astype(int)


def istriu(am):
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is not relevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return sp.all(am.row < am.col)


def istril(am):
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is not relevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return sp.all(am.row > am.col)


def istriangular(am):
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return istril(am) or istriu(am)


def issymmetric(am):
    if am.shape[0] != am.shape[1]:
        logger.warning('Matrix is not square, symmetrical is not relevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    if istril(am) or istriu(am):
        return False
    inds_up = am.row < am.col
    inds_lo = am.row > am.col
    # Check if locations of non-zeros are symmetrically distributed
    if inds_up.size != inds_lo.size:
        return False
    # Then check if values are also symmetric
    data_up = am.data[inds_up]
    data_lo = am.data[inds_lo]
    return sp.allclose(data_up, data_lo)


def am_to_im(am):
    if am.shape[0] != am.shape[1]:
        raise Exception('Adjacency matrices must be square')
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    conn = sp.vstack((am.row, am.col)).T
    row = conn[:, 0]
    data = am.data
    col = sp.arange(sp.size(am.data))
    if istriangular(am):
        row = sp.append(row, conn[:, 1])
        data = sp.append(data, data)
        col = sp.append(col, col)
    im = sprs.coo.coo_matrix((data, (row, col)), (row.max()+1, col.max()+1))
    return im


def im_to_am(im):
    if im.shape[0] == im.shape[1]:
        print('Warning: Received matrix is square which is unlikely')
    if im.shape[0] > im.shape[1]:
        print('Warning: Received matrix has more sites than bonds')
    if im.format != 'coo':
        im = im.tocoo(copy=False)


def trim(network, pores=[], throats=[]):
    '''
    Remove pores or throats from the network.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network from which pores or throats should be removed

    pores (or throats) : array_like
        The indices of the of the pores or throats to be removed from the
        network.

    Notes
    -----
    This is an in-place operation, meaning the received Network object will
    be altered directly.


    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn.Np
    125
    >>> pn.Nt
    300
    >>> op.topotools.trim(network=pn, pores=[1])
    >>> pn.Np
    124
    >>> pn.Nt
    296

    '''
    pores = sp.array(pores, ndmin=1)
    throats = sp.array(throats, ndmin=1)
    Pkeep = sp.copy(network['pore.all'])
    Tkeep = sp.copy(network['throat.all'])
    if sp.size(pores) > 0:
        Pkeep[pores] = False
        if not sp.any(Pkeep):
            raise Exception('Cannot delete ALL pores')
        # Performing customized find_neighbor_throats which is much faster, but
        # not general for other types of queries
        temp = sp.in1d(network['throat.conns'].flatten(), pores)
        temp = sp.reshape(temp, (network.Nt, 2))
        Ts = sp.any(temp, axis=1)
        if len(Ts) > 0:
            Tkeep[Ts] = False
    if sp.size(throats) > 0:
        Tkeep[throats] = False
        # The folling IF catches the special case of deleting ALL throats
        # It removes all throat props, adds 'all', and skips rest of function
        if not sp.any(Tkeep):
            logger.info('Removing ALL throats from network')
            for item in network.keys():
                if item.split('.')[0] == 'throat':
                    del network[item]
            network['throat.all'] = sp.array([], ndmin=1)
            return

    # Temporarily store throat conns and pore map for processing later
    Np_old = network.Np
    Nt_old = network.Nt
    Pkeep_inds = sp.where(Pkeep)[0]
    Tkeep_inds = sp.where(Tkeep)[0]
    Pmap = sp.ones((network.Np,), dtype=int)*-1
    Pid = network['pore._id']
    Tid = network['throat._id']
    tpore1 = network['throat.conns'][:, 0]
    tpore2 = network['throat.conns'][:, 1]

    # Delete specified pores and throats from all objects
    for obj in network.project:
        if (obj.Np == Np_old) and (obj.Nt == Nt_old):
            Ps = Pkeep_inds
            Ts = Tkeep_inds
        else:
            Ps = obj.map_pores(ids=Pid[Pkeep])
            Ts = obj.map_throats(ids=Tid[Tkeep])
        for key in list(obj.keys()):
            temp = obj.pop(key)
            if key.split('.')[0] == 'throat':
                obj.update({key: temp[Ts]})
            if key.split('.')[0] == 'pore':
                obj.update({key: temp[Ps]})

    # Remap throat connections
    Pmap[Pkeep] = sp.arange(0, sp.sum(Pkeep))
    Tnew1 = Pmap[tpore1[Tkeep]]
    Tnew2 = Pmap[tpore2[Tkeep]]
    network.update({'throat.conns': sp.vstack((Tnew1, Tnew2)).T})

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def extend(network, pore_coords=[], throat_conns=[], labels=[]):
    r'''
    Add individual pores and/or throats to the network from a list of coords
    or conns.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network to which pores or throats should be added

    pore_coords : array_like
        The coordinates of the pores to add

    throat_conns : array_like
        The throat connections to add

    labels : string, or list of strings, optional
        A list of labels to apply to the new pores and throats

    Notes
    -----
    This needs to be enhanced so that it increases the size of all pore
    and throat props and labels on ALL associated Phase objects.  At the
    moment it throws an error is there are any associated Phases.

    This is an in-place operation, meaning the received Network object will
    be altered directly.

    '''
    if len(network.project.phases()) > 0:
        raise Exception('Project has active Phases, copy network to a new ' +
                        'project and try again')

    Np_old = network.num_pores()
    Nt_old = network.num_throats()
    Np = Np_old + int(sp.size(pore_coords)/3)
    Nt = Nt_old + int(sp.size(throat_conns)/2)
    # Adjust 'all' labels
    del network['pore.all'], network['throat.all']
    network['pore.all'] = sp.ones((Np,), dtype=bool)
    network['throat.all'] = sp.ones((Nt,), dtype=bool)
    # Add coords and conns
    if sp.size(pore_coords) > 0:
        coords = sp.vstack((network['pore.coords'], pore_coords))
        network['pore.coords'] = coords
    if sp.size(throat_conns) > 0:
        conns = sp.vstack((network['throat.conns'], throat_conns))
        network['throat.conns'] = conns

    # Increase size of any prop or label arrays aready on Network
    for item in list(network.keys()):
        if item.split('.')[1] not in ['coords', 'conns', 'all', '_id']:
            N = network._count(element=item.split('.')[0])
            arr = network.pop(item)
            s = arr.shape
            network[item] = sp.zeros(shape=(N, *s[1:]), dtype=arr.dtype)
            # This is a temporary work-around until I learn to handle 2+ dims
            network[item][:arr.shape[0]] = arr

    # Apply labels, if supplied
    if labels != []:
        # Convert labels to list if necessary
        if type(labels) is str:
            labels = [labels]
        for label in labels:
            # Remove pore or throat from label, if present
            label = label.split('.')[-1]
            if sp.size(pore_coords) > 0:
                Ps = sp.r_[Np_old:Np]
                if 'pore.'+label not in network.labels():
                    network['pore.'+label] = False
                network['pore.'+label][Ps] = True
            if sp.size(throat_conns) > 0:
                Ts = sp.r_[Nt_old:Nt]
                if 'throat.'+label not in network.labels():
                    network['throat.'+label] = False
                network['throat.'+label][Ts] = True


def reduce_coordination(network, z):
    r"""
    """
    # Find minimum spanning tree using random weights
    am = network.create_adjacency_matrix(weights=sp.rand(network.Nt),
                                         triu=False)
    mst = csgraph.minimum_spanning_tree(am, overwrite=True)
    mst = mst.tocoo()

    # Label throats on spanning tree to avoid deleting them
    Ts = network.find_connecting_throat(mst.row, mst.col)
    Ts = sp.hstack(Ts)
    network['throat.mst'] = False
    network['throat.mst'][Ts] = True

    # Trim throats not on the spanning tree to acheive desired coordination
    Ts = sp.random.permutation(network.throats('mst', mode='complement'))
    Ts = Ts[:int(network.Nt - network.Np*(z/2))]
    trim(network=network, throats=Ts)


def label_faces(network, tol=0.1):
    r"""
    Finds pores on the surface of the network and labels them according to
    whether they are on the *top*, *bottom*, etc.  This function assumes the
    network is cubic in shape (i.e. with six flat sides)

    Parameters
    ----------
    network : OpenPNM Network object
        The network to apply the labels
    """
    find_surface_pores(network)
    Psurf = network['pore.surface']
    crds = network['pore.coords']
    xmin, xmax = sp.amin(crds[:, 0]), sp.amax(crds[:, 0])
    xspan = xmax - xmin
    ymin, ymax = sp.amin(crds[:, 1]), sp.amax(crds[:, 1])
    yspan = ymax - ymin
    zmin, zmax = sp.amin(crds[:, 2]), sp.amax(crds[:, 2])
    zspan = zmax - zmin
    network['pore.back'] = (crds[:, 0] > (1-tol)*xmax) * Psurf
    network['pore.right'] = (crds[:, 1] > (1-tol)*ymax) * Psurf
    network['pore.top'] = (crds[:, 2] > (1-tol)*zmax) * Psurf
    network['pore.front'] = (crds[:, 0] < (xmin + tol*xspan)) * Psurf
    network['pore.left'] = (crds[:, 1] < (xmin + tol*xspan)) * Psurf
    network['pore.bottom'] = (crds[:, 2] < (xmin + tol*xspan)) * Psurf


def find_surface_pores(network, markers=None, label='surface'):
    r"""
    Find the pores on the surface of the domain by performing a Delaunay
    triangulation between the network pores and some external ``markers``. All
    pores connected to these external marker points are considered surface
    pores.

    Parameters
    ----------
    network: OpenPNM Network Object
        The network for which the surface pores are to be found

    markers: array_like
        3 x N array of the marker coordiantes to use in the triangulation.  The
        labeling is performed in one step, so all points are added, and then
        any pores connected to at least one marker is given the provided label.
        By default, this function will automatically generate 6 points outside
        each axis of the network domain.

        Users may wish to specify a single external marker point and provide an
        appropriate label in order to identify specific faces.  For instance,
        the marker may be *above* the domain, and the label might be
        'top_surface'.

    label : string
        The label to apply to the pores.  The default is 'surface'.

    Notes
    -----
    This function does not check whether the given markers actually lie outside
    the domain, allowing the labeling of *internal* sufaces.

    If this method fails to mark some surface pores, consider sending more
    markers on each face.

    Examples
    --------
    >>> import openpnm as op
    >>> net = op.network.Cubic(shape=[5, 5, 5])
    >>> op.topotools.find_surface_pores(network=net)
    >>> net.num_pores('surface')
    98

    When cubic networks are created, the surfaces are already labeled:

    >>> net.num_pores(['top','bottom', 'left', 'right', 'front','back'])
    98

    This function is mostly useful for unique networks such as spheres, random
    topology, or networks that have been subdivied.

    """
    import scipy.spatial as sptl
    if markers is None:
        (xmax, ymax, zmax) = sp.amax(network['pore.coords'], axis=0)
        (xmin, ymin, zmin) = sp.amin(network['pore.coords'], axis=0)
        xave = (xmin+xmax)/2
        yave = (ymin+ymax)/2
        zave = (zmin+zmax)/2
        markers = [[xmax + xave, yave, zave],
                   [xmin - xave, yave, zave],
                   [xave, ymax + yave, zave],
                   [xave, ymin - yave, zave],
                   [xave, yave, zmax + zave],
                   [xave, yave, zmin - zave]]
    markers = sp.atleast_2d(markers)
    tri = sptl.Delaunay(network['pore.coords'], incremental=True)
    tri.add_points(markers)
    (indices, indptr) = tri.vertex_neighbor_vertices
    for k in range(network.Np, tri.npoints):
        neighbors = indptr[indices[k]:indices[k+1]]
        inds = sp.where(neighbors < network.Np)
        neighbors = neighbors[inds]
        if 'pore.'+label not in network.keys():
            network['pore.'+label] = False
        network['pore.'+label][neighbors] = True


def clone_pores(network, pores, labels=['clone'], mode='parents'):
    r'''
    Clones the specified pores and adds them to the network

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which the new pores are to be added

    pores : array_like
        List of pores to clone

    labels : string, or list of strings
        The labels to apply to the clones, default is 'clone'

    mode : string
        Controls the connections between parents and clones.  Options are:

        - 'parents': (Default) Each clone is connected only to its parent
        - 'siblings': Clones are only connected to each other in the same
                      manner as parents were connected
        - 'isolated': No connections between parents or siblings
    '''
    if len(network.project.geometries()) > 0:
        logger.warning('Network has active Geometries, new pores must be \
                        assigned a Geometry')
    if len(network.project.phases()) > 0:
        raise Exception('Network has active Phases, cannot proceed')

    if type(labels) == str:
        labels = [labels]
    Np = network.Np
    Nt = network.Nt
    # Clone pores
    parents = sp.array(pores, ndmin=1)
    pcurrent = network['pore.coords']
    pclone = pcurrent[pores, :]
    pnew = sp.concatenate((pcurrent, pclone), axis=0)
    Npnew = sp.shape(pnew)[0]
    clones = sp.arange(Np, Npnew)
    # Add clone labels to network
    for item in labels:
        network['pore.'+item] = False
        network['throat.'+item] = False
    # Add connections between parents and clones
    if mode == 'parents':
        tclone = sp.vstack((parents, clones)).T
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'siblings':
        ts = network.find_neighbor_throats(pores=pores, mode='intersection')
        tclone = network['throat.conns'][ts] + network.num_pores()
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'isolated':
        extend(network=network, pore_coords=pclone)
    # Apply provided labels to cloned pores
    for item in labels:
        network['pore.'+item][network.pores('all') >= Np] = True
        network['throat.'+item][network.throats('all') >= Nt] = True

    # Any existing adjacency and incidence matrices will be invalid
    network._am.clear()
    network._im.clear()


def stitch(network, donor, P_network, P_donor, method='nearest',
           len_max=sp.inf, label_suffix=''):
    r'''
    Stitches a second a network to the current network.

    Parameters
    ----------
    networK : OpenPNM Network Object
        The Network that will to which to donor Network will be attached

    donor : OpenPNM Network Object
        The Network to stitch on to the current Network

    P_network : array_like
        The pores on the current Network

    P_donor : array_like
        The pores on the donor Network

    label_suffix : string or None
        Some text to append to each label in the donor Network before
        inserting them into the recipient.  The default is to append no
        text, but a common option would be to append the donor Network's
        name. To insert none of the donor labels, use None.

    len_max : float
        Set a length limit on length of new throats

    method : string (default = 'delaunay')
        The method to use when making pore to pore connections. Options are:

        - 'delaunay' : Use a Delaunay tessellation
        - 'nearest' : Connects each pore on the receptor network to its nearest
                      pore on the donor network

    Notes
    -----
    Before stitching it is necessary to translate the pore coordinates of
    one of the Networks so that it is positioned correctly relative to the
    other.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn2 = op.network.Cubic(shape=[5, 5, 5])
    >>> [pn.Np, pn.Nt]
    [125, 300]
    >>> [pn2.Np, pn2.Nt]
    [125, 300]
    >>> pn2['pore.coords'][:, 2] += 5.0
    >>> op.topotools.stitch(network=pn, donor=pn2, P_network=pn.pores('top'),
    ...                     P_donor=pn2.pores('bottom'), method='nearest',
    ...                     len_max=1.0)
    >>> [pn.Np, pn.Nt]
    [250, 625]

    '''
    # Ensure Networks have no associated objects yet
    if (len(network.project) > 1) or (len(donor.project) > 1):
        raise Exception('Cannot stitch a Network with active sibling objects')
    network['throat.stitched'] = False
    # Get the initial number of pores and throats
    N_init = {}
    N_init['pore'] = network.Np
    N_init['throat'] = network.Nt
    if method == 'nearest':
        P1 = P_network
        P2 = P_donor + N_init['pore']  # Increment pores on donor
        C1 = network['pore.coords'][P_network]
        C2 = donor['pore.coords'][P_donor]
        D = sp.spatial.distance.cdist(C1, C2)
        [P1_ind, P2_ind] = sp.where(D <= len_max)
        conns = sp.vstack((P1[P1_ind], P2[P2_ind])).T
    else:
        raise RuntimeError('<{}> method not supported'.format(method))

    # Enter donor's pores into the Network
    extend(network=network, pore_coords=donor['pore.coords'])

    # Enter donor's throats into the Network
    extend(network=network, throat_conns=donor['throat.conns'] +
           N_init['pore'])

    # Trim throats that are longer then given len_max
    C1 = network['pore.coords'][conns[:, 0]]
    C2 = network['pore.coords'][conns[:, 1]]
    L = sp.sum((C1 - C2)**2, axis=1)**0.5
    conns = conns[L <= len_max]

    # Add donor labels to recipient network
    if label_suffix is not None:
        if label_suffix != '':
            label_suffix = '_'+label_suffix
        for label in donor.labels():
            element = label.split('.')[0]
            locations = sp.where(network._get_indices(element) >=
                                 N_init[element])[0]
            if label + label_suffix not in network.keys():
                network[label + label_suffix] = False
            network[label+label_suffix][locations] = donor[label]

    # Add the new stitch throats to the Network
    extend(network=network, throat_conns=conns, labels='stitched')

    # Remove donor from Workspace, if present
    # This check allows for the reuse of a donor Network multiple times
    for sim in list(ws.values()):
        if donor in sim:
            del ws[sim.name]


def connect_pores(network, pores1, pores2, labels=[], add_conns=True):
    r'''
    Returns the possible connections between two group of pores, and optionally
    makes the connections.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores1 : array_like
        The first group of pores on the network

    pores2 : array_like
        The second group of pores on the network

    labels : list of strings
        The labels to apply to the new throats.  This argument is only needed
        if ``add_conns`` is True.

    add_conns : bool
        Indicates whether the connections should be added to the supplied
        network (default is True).  Otherwise, the connections are returned
        as an Nt x 2 array that can be passed directly to ``extend``.

    Notes
    -----
    It creates the connections in a format which is acceptable by
    the default OpenPNM connection ('throat.conns') and either adds them to
    the network or returns them.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn.Nt
    300
    >>> op.topotools.connect_pores(network=pn, pores1=[22, 32],
    ...                            pores2=[16, 80, 68])
    >>> pn.Nt
    306
    >>> pn['throat.conns'][300:306]
    array([[16, 22],
           [22, 80],
           [22, 68],
           [16, 32],
           [32, 80],
           [32, 68]])

    '''
    size1 = sp.size(pores1)
    size2 = sp.size(pores2)
    array1 = sp.repeat(pores1, size2)
    array2 = sp.tile(pores2, size1)
    conns = sp.vstack([array1, array2]).T
    if add_conns:
        extend(network=network, throat_conns=conns, labels=labels)
    else:
        return conns


def find_pore_to_pore_distance(network, pores1=None, pores2=None):
    r'''
    Find the distance between all pores on set one to each pore in set 2

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object containing the pore coordinates

    pores1 : array_like
        The pore indices of the first set

    pores2 : array_Like
        The pore indices of the second set.  It's OK if these indices are
        partially or completely duplicating ``pores``.

    Returns
    -------
    A distance matrix with ``len(pores1)`` rows and ``len(pores2)`` columns.
    The distance between pore *i* in ``pores1`` and *j* in ``pores2`` is
    located at *(i, j)* and *(j, i)* in the distance matrix.

    '''
    from scipy.spatial.distance import cdist
    p1 = sp.array(pores1, ndmin=1)
    p2 = sp.array(pores2, ndmin=1)
    coords = network['pore.coords']
    return cdist(coords[p1], coords[p2])


def subdivide(network, pores, shape, labels=[]):
    r'''
    It trim the pores and replace them by cubic networks with the sent shape.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores : array_like
        The first group of pores to be replaced

    shape : array_like
        The shape of cubic networks in the target locations

    Notes
    -----
    - It works only for cubic networks.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 6, 5], spacing=0.001)
    >>> pn.Np
    150
    >>> nano_pores = [2, 13, 14, 15]
    >>> op.topotools.subdivide(network=pn, pores=nano_pores, shape=[4, 7, 3],
    ...                        labels='nano')
    >>> pn.Np
    482

    '''
    mro = network._mro()
    if 'Cubic' not in mro:
        raise Exception('Subdivide is only supported for Cubic Networks')
    from openpnm.network import Cubic
    pores = sp.array(pores, ndmin=1)

    # Checks to find boundary pores in the selected pores
    if 'pore.boundary' in network.labels():
        if (sp.in1d(pores, network.pores('boundary'))).any():
            raise Exception('boundary pores cannot be subdivided!')
    if not hasattr(network, '_subdivide_flag'):
        network._subdivide_flag = True
    else:
        raise Exception('The network has subdivided pores, so the method \
                         does not support another subdivision')
    # Assigning right shape and division
    if sp.size(shape) != 2 and sp.size(shape) != 3:
        raise Exception('Subdivide not implemented for Networks other than 2D \
                         and 3D')
    elif sp.size(shape) == 3 and 1 not in shape:
        div = sp.array(shape, ndmin=1)
        single_dim = None
    else:
        single_dim = sp.where(sp.array(network._shape) == 1)[0]
        if sp.size(single_dim) == 0:
            single_dim = None
        if sp.size(shape) == 3:
            div = sp.array(shape, ndmin=1)
        else:
            div = sp.zeros(3, dtype=sp.int32)
            if single_dim is None:
                dim = 2
            else:
                dim = single_dim
            div[dim] = 1
            div[-sp.array(div, ndmin=1, dtype=bool)] = sp.array(shape,
                                                                ndmin=1)

    # Creating small network and handling labels
    networkspacing = network.spacing
    new_netspacing = networkspacing/div
    new_net = Cubic(shape=div, spacing=new_netspacing)
    main_labels = ['left', 'right', 'front', 'back', 'top', 'bottom']
    if single_dim is not None:
        label_groups = sp.array([['front', 'back'],
                                 ['left', 'right'],
                                 ['top', 'bottom']])
        non_single_labels = label_groups[sp.array([0, 1, 2]) != single_dim]
    for l in main_labels:
        new_net['pore.surface_' + l] = False
        network['pore.surface_' + l] = False
        if single_dim is None:
            new_net['pore.surface_' + l][new_net.pores(labels=l)] = True
        else:
            for ind in [0, 1]:
                loc = (non_single_labels[ind] == l)
                temp_pores = new_net.pores(non_single_labels[ind][loc])
                new_net['pore.surface_' + l][temp_pores] = True

    old_coords = sp.copy(new_net['pore.coords'])
    if labels == []:
        labels = ['pore.subdivided_' + new_net.name]
    for P in pores:
        # Shifting the new network to the right location and attaching it to
        # the main network
        shift = network['pore.coords'][P] - networkspacing/2
        new_net['pore.coords'] += shift
        Pn = network.find_neighbor_pores(pores=P)
        try:
            Pn_new_net = network.pores(labels)
        except:
            Pn_new_net = []
        Pn_old_net = Pn[~sp.in1d(Pn, Pn_new_net)]
        Np1 = network.Np
        extend(pore_coords=new_net['pore.coords'],
               throat_conns=new_net['throat.conns'] + Np1,
               labels=labels, network=network)

        # Moving the temporary labels to the big network
        for l in main_labels:
            network['pore.surface_'+l][Np1:] = new_net['pore.surface_'+l]

        # Stitching the old pores of the main network to the new extended pores
        surf_pores = network.pores('surface_*')
        surf_coord = network['pore.coords'][surf_pores]
        for neighbor in Pn:
            neighbor_coord = network['pore.coords'][neighbor]
            dist = [round(sp.inner(neighbor_coord-x, neighbor_coord-x),
                          20) for x in surf_coord]
            nearest_neighbor = surf_pores[dist == sp.amin(dist)]
            if neighbor in Pn_old_net:
                coplanar_labels = network.labels(pores=nearest_neighbor)
                new_neighbors = network.pores(coplanar_labels,
                                              mode='intersection')
                # This might happen to the edge of the small network
                if sp.size(new_neighbors) == 0:
                    labels = network.labels(pores=nearest_neighbor,
                                            mode='intersection')
                    common_label = [l for l in labels if 'surface_' in l]
                    new_neighbors = network.pores(common_label)
            elif neighbor in Pn_new_net:
                new_neighbors = nearest_neighbor
            connect_pores(network=network, pores1=neighbor,
                          pores2=new_neighbors, labels=labels)

        # Removing temporary labels
        for l in main_labels:
            network['pore.surface_' + l] = False
        new_net['pore.coords'] = sp.copy(old_coords)

    network._label_surfaces()
    for l in main_labels:
        del network['pore.surface_'+l]
    trim(network=network, pores=pores)
    ws = network.project.workspace
    ws.close_project(new_net.project)


def trim_occluded_throats(network, mask='all'):
    r"""
    Remove throats with zero area from the network and also remove
    pores that are isolated (as a result or otherwise)

    Parameters
    ----------
    network : OpenPNM Network Object

    mask : string
        Applies routine only to pores and throats with this label
    """
    occluded_ts = network['throat.area'] == 0
    if sp.sum(occluded_ts) > 0:
        occluded_ts *= network["throat."+mask]
        trim(network=network, throats=occluded_ts)


def merge_pores(network, pores, labels=['merged']):
    r"""
    Combines a selection of pores into a new single pore located at the
    centroid of the selected pores and connected to all of their neighbors.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores : array_like
        The list of pores which are to be combined into a new single pore

    labels : string or list of strings
        The labels to apply to the new pore and new throat connections

    Notes
    -----
    The selection of pores should be chosen carefully, preferrable so that
    they all form a continuous cluster.  For instance, it is recommended
    to use the ``find_nearby_pores`` method to find all pores within a
    certain distance of a given pore, and these can then be merged without
    causing any abnormal connections.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[20, 20, 1])
    >>> Ps = pn.find_nearby_pores(pores=111, r=5, flatten=True)
    >>> op.topotools.merge_pores(network=pn, pores=Ps, labels=['merged'])
    >>> print(pn.Np)
    321
    >>> pn.pores('merged')
    array([320])
    >>> pn.num_throats('merged')
    32

    """
    Pn = network.find_neighbor_pores(pores=pores,
                                     mode='union',
                                     flatten=True,
                                     excl_self=True)
    xyz = sp.mean(network['pore.coords'][pores], axis=0)
    extend(network, pore_coords=xyz, labels=labels)
    Pnew = network.Ps[-1]
    connect_pores(network, pores1=Pnew, pores2=Pn, labels=labels)
    trim(network=network, pores=pores)


def _template_sphere_disc(dim, outer_radius, inner_radius):
    r"""
    This private method generates an image array of a sphere/shell-disc/ring.
    It is useful for passing to Cubic networks as a ``template`` to make
    networks with desired shapes.

    Parameters
    ----------
    dim : int
        Network dimension

    outer_radius : int
        Number of the nodes in the outer radius of the network

    inner_radius : int
        Number of the nodes in the inner radius of the network

    Returns
    -------
    A Numpy array containing 1's to demarcate the desired shape, and 0's
    elsewhere.

    """
    rmax = sp.array(outer_radius, ndmin=1)
    rmin = sp.array(inner_radius, ndmin=1)
    ind = 2 * rmax - 1
    coord = sp.indices((ind * sp.ones(dim, dtype=int)))
    coord = coord - (ind - 1)/2
    x = coord[0, :]
    y = coord[1, :]
    if dim == 2:
        img = (x ** 2 + y ** 2) < rmax ** 2
    elif dim == 3:
        z = coord[2, :]
        img = (x ** 2 + y ** 2 + z ** 2) < rmax ** 2
    if rmin[0] != 0:
        if dim == 2:
            img_min = (x ** 2 + y ** 2) > rmin ** 2
        elif dim == 3:
            img_min = (x ** 2 + y ** 2 + z ** 2) > rmin ** 2
        img = img * img_min
    return img


def template_sphere_shell(outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a sphere-shell. It is useful for
    passing to Cubic networks as a ``template`` to make spherical shaped
    networks.

    Parameters
    ----------
    outer_radius : int
        Number of nodes in the outer radius of the sphere.

    inner_radius : int
        Number of nodes in the inner radius of the shell.  a value of 0 will
        result in a solid sphere.

    Returns
    -------
    A Numpy array containing 1's to demarcate the sphere-shell, and 0's
    elsewhere.

    """
    img = _template_sphere_disc(dim=3, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return img


def template_cylinder_annulus(height, outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a disc-ring.  It is useful for
    passing to Cubic networks as a ``template`` to make circular-shaped 2D
    networks.

    Parameters
    ----------
    height : int
        The height of the cylinder

    outer_radius : int
        Number of nodes in the outer radius of the cylinder

    inner_radius : int
        Number of the nodes in the inner radius of the annulus.  A value of 0
        will result in a solid cylinder.

    Returns
    -------
    A Numpy array containing 1's to demarcate the disc-ring, and 0's
    elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    img = sp.tile(sp.atleast_3d(img), reps=height)
    return img


def plot_connections(network, throats=None, fig=None, **kwargs):
    r"""
    Produces a 3D plot of the network topology showing how throats connect
    for quick visualization without having to export data to veiw in Paraview.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network whose topological connections to plot

    throats : array_like (optional)
        The list of throats to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no throats are
        specified then all throats are shown.

    fig and **kwargs: Matplotlib figure handle and line property arguments
        If a ``fig`` is supplied, then the topology will be overlaid.  By also
        passing in different line properties such as ``color`` and limiting
        which ``throats`` are plots, this makes it possible to plot different
        types of throats on the same plot.

        For information on available line style options, visit the Matplotlib
        documentation at:

        http://matplotlib.org/api/lines_api.html#matplotlib.lines.Line2D

    Notes
    -----
    The figure handle returned by this method can be passed into
    ``plot_coordinates`` to create a plot that combines pore coordinates and
    throat connections, and vice versa.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundary_pores()
    >>> Ts = pn.throats('*boundary', mode='complement')
    >>> # Create figure showing boundary throats
    >>> fig = op.topotools.plot_connections(network=pn, throats=Ts)
    >>> Ts = pn.throats('*boundary')
    >>> # Pass existing fig back into function to plot additional throats
    >>> fig = op.topotools.plot_connections(network=pn, throats=Ts,
    ...                                     fig=fig, color='r')

    """
    import matplotlib.pyplot as plt

    if throats is None:
        Ts = network.Ts
    else:
        Ts = network._parse_indices(indices=throats)

    ThreeD = False
    if len(sp.unique(network['pore.coords'][:, 2])) > 1:
        from mpl_toolkits.mplot3d import Axes3D
        ThreeD = True

    if fig is None:
        fig = plt.figure()
        if ThreeD:
            ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.get_axes()[0]

    # Create dummy indexing to sp.inf
    i = -1*sp.ones((sp.size(Ts)*3, ), dtype=int)
    i[0::3] = network['throat.conns'][Ts, 0]
    i[1::3] = network['throat.conns'][Ts, 1]

    # Collect coordinates and scale axes to fit
    Ps = sp.unique(network['throat.conns'][Ts])
    X = network['pore.coords'][Ps, 0]
    Y = network['pore.coords'][Ps, 1]
    Z = network['pore.coords'][Ps, 2]
    if ThreeD:
        _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)

    # Add sp.inf to the last element of pore.coords (i.e. -1)
    inf = sp.array((sp.inf,))
    X = sp.hstack([network['pore.coords'][:, 0], inf])
    Y = sp.hstack([network['pore.coords'][:, 1], inf])
    Z = sp.hstack([network['pore.coords'][:, 2], inf])
    if ThreeD:
        ax.plot(xs=X[i], ys=Y[i], zs=Z[i], **kwargs)
    else:
        plt.plot(X[i], Y[i], **kwargs)

    return fig


def plot_coordinates(network, pores=None, fig=None, **kwargs):
    r"""
    Produces a 3D plot showing specified pore coordinates as markers

    Parameters
    ----------
    network : OpenPNM Network Object
        The network whose topological connections to plot

    pores : array_like (optional)
        The list of pores to plot if only a sub-sample is desired.  This is
        useful for inspecting a small region of the network.  If no pores are
        specified then all are shown.

    fig and **kwargs: Matplotlib figure handle and line property arguments
        If a ``fig`` is supplied, then the topology will be overlaid.  By also
        passing in different marker properties such as size (``s``) and
        limiting which ``pores`` are plotted, this makes it possible to plot
        different types of pores on the same plot.

        For information on available marker style options, visit the Matplotlib
        documentation at:

        http://matplotlib.org/api/lines_api.html#matplotlib.lines.Line2D

    Notes
    -----
    The figure handle returned by this method can be passed into
    ``plot_topology`` to create a plot that combines pore coordinates and
    throat connections, and vice versa.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundary_pores()
    >>> Ps = pn.pores('internal')
    >>> # Create figure showing internal pores
    >>> fig = op.topotools.plot_coordinates(network=pn, pores=Ps, c='b')
    >>> Ps = pn.pores('*boundary')
    >>> # Pass existing fig back into function to plot boundary pores
    >>> fig = op.topotools.plot_coordinates(network=pn, pores=Ps, fig=fig,
    ...                                         c='r')

    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if pores is None:
        Ps = network.Ps
    else:
        Ps = network._parse_indices(indices=pores)

    ThreeD = False
    if len(sp.unique(network['pore.coords'][:, 2])) > 1:
        from mpl_toolkits.mplot3d import Axes3D
        ThreeD = True

    if fig is None:
        fig = plt.figure()
        if ThreeD:
            ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.get_axes()[0]

    # Collect specified coordinates
    X = network['pore.coords'][Ps, 0]
    Y = network['pore.coords'][Ps, 1]
    Z = network['pore.coords'][Ps, 2]
    _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)

    if ThreeD:
        ax.scatter(xs=X, ys=Y, zs=Z, **kwargs)
    else:
        plt.scatter(X, Y, **kwargs)
    return fig


def _scale_3d_axes(ax, X, Y, Z):
    if hasattr(ax, '_scaled'):
        logger.warning('Axes is already scaled to previously plotted data')
    else:
        ax._scaled = True
        max_range = sp.array([X.max()-X.min(), Y.max()-Y.min(),
                              Z.max()-Z.min()]).max() / 2.0
        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


def generate_base_points(num_points, domain_size, prob=None):
    r"""
    Generates a set of base points for passing into the DelaunayVoronoiDual
    class.  The points can be distributed in spherical, cylindrical, or
    rectilinear patterns.

    Parameters
    ----------
    num_points : scalar
        The number of base points that lie within the domain.  Note that the
        actual number of points returned will be larger, with the extra points
        lying outside the domain.

    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].

    prob : 3D array, optional
        A 3D array that contains fractional (0-1) values indicating the
        liklihood that a point in that region should be kept.  If not specified
        an array containing 1's in the shape of a sphere, cylinder, or cube is
        generated, depnending on the give ``domain_size`` with zeros outside.
        When specifying a custom probabiliy map is it recommended to also set
        values outside the given domain to zero.  If not, then the correct
        shape will still be returned, but with too few points in it.

    Notes
    -----
    This method places the given number of points within the specified domain,
    then reflects these points across each domain boundary.  This results in
    smooth flat faces at the boundaries once these excess pores are trimmed.

    The reflection approach tends to create larger pores near the surfaces, so
    it might be necessary to use the ``prob`` argument to specify a slightly
    higher density of points near the surfaces.

    For rough faces, it is necessary to define a larger than desired domain
    then trim to the desired size.  This will discard the reflected points
    plus some of the original points.

    Examples
    --------
    The following generates a spherical array with higher values near the core.
    It uses a distance transform to create a sphere of radius 10, then a
    second distance transform to create larger values in the center away from
    the sphere surface.  These distance values could be further skewed by
    applying a power, with values higher than 1 resulting in higher values in
    the core, and fractional values smoothinging them out a bit.

    >>> import openpnm as op
    >>> import scipy as sp
    >>> import scipy.ndimage as spim
    >>> im = sp.ones([21, 21, 21], dtype=int)
    >>> im[10, 10, 10] = 0
    >>> im = spim.distance_transform_edt(im) <= 20  # Create sphere of 1's
    >>> prob = spim.distance_transform_edt(im)
    >>> prob = prob / sp.amax(prob)  # Normalize between 0 and 1
    >>> pts = op.topotools.generate_base_points(num_points=50,
    ...                                         domain_size=[2],
    ...                                         prob=prob)
    >>> net = op.network.DelaunayVoronoiDual(points=pts, shape=[2])
    """
    def _try_points(num_points, prob):
        prob = sp.array(prob)/sp.amax(prob)  # Ensure prob is normalized
        base_pts = []
        N = 0
        while N < num_points:
            pt = sp.random.rand(3)  # Generate a point
            # Test whether to keep it or not
            [indx, indy, indz] = sp.floor(pt*sp.shape(prob)).astype(int)
            if sp.random.rand(1) <= prob[indx][indy][indz]:
                base_pts.append(pt)
                N += 1
        base_pts = sp.array(base_pts)
        return base_pts
    if len(domain_size) == 1:  # Spherical
        domain_size = sp.array(domain_size)
        if prob is None:
            prob = sp.ones([41, 41, 41])
            prob[20, 20, 20] = 0
            prob = spim.distance_transform_bf(prob) <= 20
        base_pts = _try_points(num_points, prob)
        # Convert to spherical coordinates
        [X, Y, Z] = sp.array(base_pts - [0.5, 0.5, 0.5]).T  # Center at origin
        r = 2*sp.sqrt(X**2 + Y**2 + Z**2)*domain_size[0]
        theta = 2*sp.arctan(Y/X)
        phi = 2*sp.arctan(sp.sqrt(X**2 + Y**2)/Z)
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, phi] = [r[inds], theta[inds], phi[inds]]
        # Reflect base points across perimeter
        new_r = 2*domain_size - r
        r = sp.hstack([r, new_r])
        theta = sp.hstack([theta, theta])
        phi = sp.hstack([phi, phi])
        # Convert to Cartesean coordinates
        X = r*sp.cos(theta)*sp.sin(phi)
        Y = r*sp.sin(theta)*sp.sin(phi)
        Z = r*sp.cos(phi)
        base_pts = sp.vstack([X, Y, Z]).T
    elif len(domain_size) == 2:  # Cylindrical
        domain_size = sp.array(domain_size)
        if prob is None:
            prob = sp.ones([41, 41, 41])
            prob[20, 20, :] = 0
            prob = spim.distance_transform_bf(prob) <= 20
        base_pts = _try_points(num_points, prob)
        # Convert to cylindrical coordinates
        [X, Y, Z] = sp.array(base_pts - [0.5, 0.5, 0]).T  # Center on z-axis
        r = 2*sp.sqrt(X**2 + Y**2)*domain_size[0]
        theta = 2*sp.arctan(Y/X)
        z = Z*domain_size[1]
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        inds = ~((z > domain_size[1]) + (z < 0))
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        # Reflect base points about faces and perimeter
        new_r = 2*domain_size[0] - r
        r = sp.hstack([r, new_r])
        theta = sp.hstack([theta, theta])
        z = sp.hstack([z, z])
        r = sp.hstack([r, r, r])
        theta = sp.hstack([theta, theta, theta])
        z = sp.hstack([z, -z, 2-z])
        # Convert to Cartesean coordinates
        X = r*sp.cos(theta)
        Y = r*sp.sin(theta)
        Z = z
        base_pts = sp.vstack([X, Y, Z]).T
    elif len(domain_size) == 3:  # Rectilinear
        if prob is None:
            prob = sp.ones([10, 10, 10], dtype=float)
        base_pts = _try_points(num_points, prob)
        base_pts = base_pts*domain_size
        # Add reflected points
        base_pts = reflect_base_points(base_pts, domain_size)
    return base_pts


def reflect_base_points(base_pts=None, domain_size=None):
    r'''
    Helper function for relecting a set of points about the faces of a
    rectangular domain

    Parameters
    ----------
    base_pts : 3d array
        The coordinates of the base_pts to be reflected

    domain_size : list or array of length 3
        The upper coordinate of the face normal to the reflection along
        each axis. Lower bound is assumed to be zero
    '''
    domain_size = sp.array(domain_size)
    Nx, Ny, Nz = domain_size
    # Reflect base points about all 6 faces
    orig_pts = base_pts
    base_pts = sp.vstack((base_pts, [-1, 1, 1]*orig_pts +
                                    [2.0*Nx, 0, 0]))
    base_pts = sp.vstack((base_pts, [1, -1, 1]*orig_pts +
                                    [0, 2.0*Ny, 0]))
    base_pts = sp.vstack((base_pts, [1, 1, -1]*orig_pts +
                                    [0, 0, 2.0*Nz]))
    base_pts = sp.vstack((base_pts, [-1, 1, 1]*orig_pts))
    base_pts = sp.vstack((base_pts, [1, -1, 1]*orig_pts))
    base_pts = sp.vstack((base_pts, [1, 1, -1]*orig_pts))
    return base_pts


def find_clusters(network, mask=[], t_labels=False):
    r"""
    Identify connected clusters of pores in the network.  This method can
    also return a list of throat cluster numbers, which correspond to the
    cluster numbers of the pores to which the throat is connected.  Either
    site and bond percolation can be considered, see description of input
    arguments for details.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network

    mask : array_like, boolean
        A list of active bonds or sites (throats or pores).  If the mask is
        Np long, then the method will perform a site percolation, and if
        the mask is Nt long bond percolation will be performed.

    Returns
    -------
    A tuple containing an Np long list of pore cluster labels, and an Nt-long
    list of throat cluster labels.  The label numbers correspond such that
    pores and throats with the same label are part of the same cluster.

    Examples
    --------
    >>> import openpnm as op
    >>> from scipy import rand
    >>> pn = op.network.Cubic(shape=[25, 25, 1])
    >>> pn['pore.seed'] = rand(pn.Np)
    >>> pn['throat.seed'] = rand(pn.Nt)


    """
    # Parse the input arguments
    mask = sp.array(mask, ndmin=1)
    if mask.dtype != bool:
        raise Exception('Mask must be a boolean array of Np or Nt length')

    # If pore mask was given perform site percolation
    if sp.size(mask) == network.Np:
        (p_clusters, t_clusters) = _site_percolation(network, mask)
    # If pore mask was given perform bond percolation
    elif sp.size(mask) == network.Nt:
        (p_clusters, t_clusters) = _bond_percolation(network, mask)
    else:
        raise Exception('Mask received was neither Nt nor Np long')

    return (p_clusters, t_clusters)


def _site_percolation(network, pmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Find throats that produce site percolation
    conns = sp.copy(network['throat.conns'])
    conns[:, 0] = pmask[conns[:, 0]]
    conns[:, 1] = pmask[conns[:, 1]]
    # Only if both pores are True is the throat set to True
    tmask = sp.all(conns, axis=1)

    # Perform the clustering using scipy.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Adjust cluster numbers such that non-invaded pores are labelled -1
    # Note: The following line also takes care of assigning cluster numbers
    # to single isolated invaded pores
    p_clusters = (clusters + 1)*(pmask) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']]
    ind = (t_clusters[:, 0] == t_clusters[:, 1])
    t_clusters = t_clusters[:, 0]
    # Label non-invaded throats with -1
    t_clusters[~ind] = -1

    return (p_clusters, t_clusters)


def _bond_percolation(network, tmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Perform the clustering using scipy.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Convert clusters to a more usable output:
    # Find pores attached to each invaded throats
    Ps = network.find_connected_pores(throats=tmask, flatten=True)
    # Adjust cluster numbers such that non-invaded pores are labelled -1
    p_clusters = (clusters + 1)*(network.tomask(pores=Ps).astype(int)) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']][:, 0]
    # Label non-invaded throats with -1
    t_clusters[~tmask] = -1

    return (p_clusters, t_clusters)


def add_boundary_pores(network, pores, offset, apply_label='boundary'):
    r"""
    This method uses ``clone_pores`` to clone the input pores, then shifts
    them the specified amount and direction, then applies the given label.

    Parameters
    ----------
    pores : array_like
        List of pores to offset.  If no pores are specified, then it
        assumes that all surface pores are to be cloned.

    offset : 3 x 1 array
        The distance in vector form which the cloned boundary pores should
        be offset.

    apply_label : string
        This label is applied to the boundary pores.  Default is
        'boundary'.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> print(pn.Np)  # Confirm initial Network size
    125
    >>> Ps = pn.pores('top')  # Select pores on top face
    >>> pn.add_boundary_pores(labels=['top'])
    >>> print(pn.Np)  # Confirm addition of 25 new pores
    150

    """
    # Parse the input pores
    Ps = sp.array(pores, ndmin=1)
    if Ps.dtype is bool:
        Ps = network.toindices(Ps)
    if sp.size(pores) == 0:  # Handle an empty array if given
        return sp.array([], dtype=sp.int64)
    # Clone the specifed pores
    clone_pores(network=network, pores=Ps)
    newPs = network.pores('pore.clone')
    del network['pore.clone']
    # Offset the cloned pores
    network['pore.coords'][newPs] += offset
    # Apply labels to boundary pores (trim leading 'pores' if present)
    label = apply_label.split('.')[-1]
    label = 'pore.' + label
    network[label] = False
    network[label][newPs] = True


def find_path(network, pore_pairs, weights=None):
    r"""
    Find the shortest path between pairs of pores.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object on which the search should be performed

    pore_pairs : array_like
        An N x 2 array containing N pairs of pores for which the shortest
        path is sought.

    weights : array_like, optional
        An Nt-long list of throat weights for the search.  Typically this
        would be the throat lengths, but could also be used to represent
        the phase configuration.  If no weights are given then the
        standard topological connections of the Network are used.

    Returns
    -------
    A dictionary containing both the pores and throats that define the
    shortest path connecting each pair of input pores.

    Notes
    -----
    The shortest path is found using Dijkstra's algorithm included in the
    scipy.sparse.csgraph module

    TODO: The returned throat path contains the correct values, but not
    necessarily in the true order

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[3, 3, 3])
    >>> a = op.topotools.find_path(network=pn, pore_pairs=[[0, 4], [0, 10]])
    >>> a['pores']
    [array([0, 1, 4]), array([ 0,  1, 10])]
    >>> a['throats']
    [array([ 0, 19]), array([ 0, 37])]
    """
    Ps = sp.array(pore_pairs, ndmin=2)
    if weights is None:
        weights = sp.ones_like(network.Ts)
    graph = network.create_adjacency_matrix(weights=weights, fmt='csr',
                                            drop_zeros=False)
    paths = csgraph.dijkstra(csgraph=graph, indices=Ps[:, 0],
                             return_predecessors=True)[1]
    pores = []
    throats = []
    for row in range(0, sp.shape(Ps)[0]):
        j = Ps[row][1]
        ans = []
        while paths[row][j] > -9999:
            ans.append(j)
            j = paths[row][j]
        ans.append(Ps[row][0])
        ans.reverse()
        pores.append(sp.array(ans, dtype=int))
        Ts = network.find_neighbor_throats(pores=ans, mode='intersection')
        throats.append(sp.array(Ts, dtype=int))
    pdict = PrintableDict
    dict_ = pdict({'pores': pores, 'throats': throats})
    return dict_


def iscoplanar(coords):
    r'''
    Determines if given pores are coplanar with each other

    Parameters
    ----------
    coords : array_like
        List of pore coords to check for coplanarity.  At least 3 pores are
        required.

    Returns
    -------
    A boolean value of whether given points are coplanar (True) or not (False)
    '''
    coords = sp.array(coords, ndmin=1)
    if sp.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:, 0]
    Py = coords[:, 1]
    Pz = coords[:, 2]

    # Do easy check first, for common coordinate
    if sp.shape(sp.unique(Px))[0] == 1:
        return True
    if sp.shape(sp.unique(Py))[0] == 1:
        return True
    if sp.shape(sp.unique(Pz))[0] == 1:
        return True

    # Perform rigorous check using vector algebra
    n1 = sp.array((Px[1] - Px[0], Py[1] - Py[0], Pz[1] - Pz[0])).T
    n2 = sp.array((Px[2] - Px[1], Py[2] - Py[1], Pz[2] - Pz[1])).T
    n = sp.cross(n1, n2)
    r = sp.array((Px[1:-1] - Px[0], Py[1:-1] - Py[0], Pz[1:-1] - Pz[0]))

    n_dot = sp.dot(n, r)

    if sp.sum(n_dot) == 0:
        return True
    else:
        return False
