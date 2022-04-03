import numpy as np
import scipy as sp


def trim(network, pores=[], throats=[]):
    """
    Remove pores or throats from the network

    Parameters
    ----------
    network : GenericNetwork
        The Network from which pores or throats should be removed
    pores (or throats) : array_like
        The indices of the of the pores or throats to be removed from the
        network.

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

    """
    pores = network._parse_indices(pores)
    throats = network._parse_indices(throats)
    Pkeep = np.copy(network['pore.all'])
    Tkeep = np.copy(network['throat.all'])
    if np.size(pores) > 0:
        Pkeep[pores] = False
        if not np.any(Pkeep):
            raise Exception('Cannot delete ALL pores')
        # # Performing customized find_neighbor_throats which is much faster, but
        # # not general for other types of queries
        # temp = np.in1d(network['throat.conns'].flatten(), pores)
        # temp = np.reshape(temp, (network.Nt, 2))
        # Ts = np.any(temp, axis=1)
        # Ts = network.Ts[Ts]
        Ts = network.find_neighbor_throats(pores=~Pkeep, mode='union')
        if len(Ts) > 0:
            Tkeep[Ts] = False
    if np.size(throats) > 0:
        Tkeep[throats] = False
        # The following IF catches the special case of deleting ALL throats
        # It removes all throat props, adds 'all', and skips rest of function
        if not np.any(Tkeep):
            for item in list(network.keys()):
                if item.split('.')[0] == 'throat':
                    del network[item]
            network['throat.all'] = np.array([], ndmin=1)
            return

    # Temporarily store throat conns and pore map for processing later
    Np_old = network.Np
    Nt_old = network.Nt
    Pkeep_inds = np.where(Pkeep)[0]
    Tkeep_inds = np.where(Tkeep)[0]
    Pmap = np.ones((network.Np,), dtype=int)*-1
    tpore1 = network['throat.conns'][:, 0]
    tpore2 = network['throat.conns'][:, 1]

    # Delete specified pores and throats from all objects
    for obj in network.project[::-1]:
        if (obj.Np == Np_old) and (obj.Nt == Nt_old):
            Ps = Pkeep_inds
            Ts = Tkeep_inds
        else:  # If subdomain object then Np/Nt < Np/Nt_old
            Ps = obj.to_local(pores=Pkeep_inds, missing_vals=None)
            Ts = obj.to_local(throats=Tkeep_inds, missing_vals=None)
        for key in list(obj.keys()):
            temp = obj.pop(key)
            if key.split('.')[0] == 'throat':
                obj.update({key: temp[Ts]})
            if key.split('.')[0] == 'pore':
                obj.update({key: temp[Ps]})

    # Remap throat connections
    Pmap[Pkeep] = np.arange(0, np.sum(Pkeep))
    Tnew1 = Pmap[tpore1[Tkeep]]
    Tnew2 = Pmap[tpore2[Tkeep]]
    network.update({'throat.conns': np.vstack((Tnew1, Tnew2)).T})

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def extend(network, coords=[], conns=[], labels=[], **kwargs):
    r"""
    Add pores or throats to the network from a list of coords or conns.

    Parameters
    ----------
    network : GenericNetwork
        The network to which pores or throats should be added
    coords : array_like
        The coordinates of the pores to add.  These will be appended to the
        'pore.coords' array so should be of shape N-by-3, where N is the
        number of pores in the list.
    conns : array_like
        The throat connections to add.  These will be appended to the
        'throat.conns' array so should be of shape N-by-2.  Note that the
        numbering must point to existing pores.
    labels : str, or list[str], optional
        A list of labels to apply to the new pores and throats

    """
    if 'throat_conns' in kwargs.keys():
        conns = kwargs['throat_conns']
    if 'pore_coords' in kwargs.keys():
        coords = kwargs['pore_coords']
    coords = np.array(coords)
    conns = np.array(conns)
    Np_old = network.num_pores()
    Nt_old = network.num_throats()
    Np = Np_old + coords.shape[0]
    Nt = Nt_old + conns.shape[0]
    if np.any(conns > Np):
        raise Exception('Some throat conns point to non-existent pores')
    network.update({'pore.all': np.ones([Np, ], dtype=bool),
                    'throat.all': np.ones([Nt, ], dtype=bool)})
    # Add coords and conns
    if np.size(coords) > 0:
        coords = np.vstack((network['pore.coords'], coords))
        network['pore.coords'] = coords
    if np.size(conns) > 0:
        conns = np.vstack((network['throat.conns'], conns))
        network['throat.conns'] = conns

    # Increase size of any prop or label arrays already on network and phases
    objs = list(network.project.phases().values())
    objs.append(network)
    for obj in objs:
        obj.update({'pore.all': np.ones([Np, ], dtype=bool),
                    'throat.all': np.ones([Nt, ], dtype=bool)})
        for item in list(obj.keys()):
            N = obj._count(element=item.split('.')[0])
            if obj[item].shape[0] < N:
                arr = obj.pop(item)
                s = arr.shape
                if arr.dtype == bool:
                    obj[item] = np.zeros(shape=(N, *s[1:]), dtype=bool)
                else:
                    obj[item] = np.ones(shape=(N, *s[1:]), dtype=float)*np.nan
                obj[item][:arr.shape[0]] = arr

    # Regenerate models on all objects to fill new elements
    for obj in network.project.phases().values():
        if hasattr(obj, 'models'):
            obj.regenerate_models()

    # Apply labels, if supplied
    if labels != []:
        # Convert labels to list if necessary
        if isinstance(labels, str):
            labels = [labels]
        for label in labels:
            # Remove pore or throat from label, if present
            label = label.split('.')[-1]
            if np.size(coords) > 0:
                Ps = np.r_[Np_old:Np]
                if 'pore.'+label not in network.labels():
                    network['pore.'+label] = False
                network['pore.'+label][Ps] = True
            if np.size(conns) > 0:
                Ts = np.r_[Nt_old:Nt]
                if 'throat.'+label not in network.labels():
                    network['throat.'+label] = False
                network['throat.'+label][Ts] = True

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def clone_pores(network, pores, labels=['clone'], mode='parents'):
    r"""
    Clones the specified pores and adds them to the network

    Parameters
    ----------
    network : GenericNetwork
        The Network object to which the new pores are to be added
    pores : array_like
        List of pores to clone
    labels : str, or list[str]
        The labels to apply to the clones, default is 'clone'
    mode : str
        Controls the connections between parents and clones.  Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'parents'    Each clone is connected only to its parent.(Default)
            'siblings'   Clones are only connected to each other in the same
                         manner as parents were connected.
            'isolated'   No connections between parents or siblings
            ===========  =====================================================

    """
    if isinstance(labels, str):
        labels = [labels]
    network._parse_indices(pores)
    Np = network.Np
    Nt = network.Nt
    # Clone pores
    parents = np.array(pores, ndmin=1)
    pcurrent = network['pore.coords']
    pclone = pcurrent[pores, :]
    pnew = np.concatenate((pcurrent, pclone), axis=0)
    Npnew = np.shape(pnew)[0]
    clones = np.arange(Np, Npnew)
    # Create cloned pores first
    extend(network=network, pore_coords=pclone)
    # Apply provided labels to cloned pores
    for item in labels:
        network.set_label(label=item, pores=range(Np, Npnew))
    # Add connections between parents and clones
    if mode == 'parents':
        tclone = np.vstack((parents, clones)).T
        extend(network=network, conns=tclone)
    elif mode == 'siblings':
        ts = network.find_neighbor_throats(pores=pores, mode='xnor')
        mapping = np.zeros([network.Np, ], dtype=int)
        mapping[pores] = np.arange(Np, network.Np)
        tclone = mapping[network['throat.conns'][ts]]
        extend(network=network, throat_conns=tclone)
    elif mode == 'isolated':
        pass
    Ntnew = network.Nt
    for item in labels:
        network.set_label(label=item, throats=range(Nt, Ntnew))

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def merge_networks(network, donor=[]):
    r"""
    Combine multiple networks into one

    This does not attempt any topological manipulations (such as stiching
    nearby pores to each other).

    Parameters
    ----------
    network : GenericNetwork
        The network to which all the other networks should be added.
    donor : GenericNetwork or list of Objects
        The network object(s) to add to the given network

    See Also
    --------
    extend
    trim
    stitch

    """
    if isinstance(donor, list):
        donors = donor
    else:
        donors = [donor]

    # First fix up geometries
    # main_proj = network.project
    # main_geoms = main_proj.geometries()
    for donor in donors:
        proj = donor.project
        geoms = proj.geometries().values()
        for geo in geoms:
            if geo.name in network.project.names:
                geo.name = network.project._generate_name(geo)
            network.project.append(geo)

    for donor in donors:
        network['pore.coords'] = np.vstack((network['pore.coords'],
                                            donor['pore.coords']))
        network['throat.conns'] = np.vstack((network['throat.conns'],
                                             donor['throat.conns']
                                             + network.Np))
        p_all = np.ones((np.shape(network['pore.coords'])[0],), dtype=bool)
        t_all = np.ones((np.shape(network['throat.conns'])[0],), dtype=bool)
        network.update({'pore.all': p_all})
        network.update({'throat.all': t_all})
        for key in set(network.keys()).union(set(donor.keys())):
            if key.split('.')[1] not in ['conns', 'coords', '_id', 'all']:
                if key in network.keys():
                    pop_flag = False
                    # If key not on donor add it first with dummy values to
                    # simplify merging later
                    if key not in donor.keys():
                        if network[key].dtype == bool:  # Deal with labels
                            donor[key] = False
                        else:  # Deal with numerical data
                            element = key.split('.')[0]
                            shape = list(network[key].shape)
                            N = donor._count(element)
                            shape[0] = N
                            donor[key] = np.empty(shape=shape)*np.nan
                        pop_flag = True
                    # Then merge it with existing array on network
                    if len(network[key].shape) == 1:
                        temp = np.hstack((network[key], donor[key]))
                    else:
                        temp = np.vstack((network[key], donor[key]))
                    network[key] = temp
                    if pop_flag:
                        donor.pop(key, None)
                else:
                    # If key not on network add it first
                    if donor[key].dtype == bool:
                        network[key] = False
                    else:
                        data_shape = list(donor[key].shape)
                        pore_prop = True if key.split(".")[0] == "pore" else False
                        data_shape[0] = network.Np if pore_prop else network.Nt
                        network[key] = np.empty(data_shape) * np.nan
                    # Then append donor values to network
                    s = np.shape(donor[key])[0]
                    network[key][-s:] = donor[key]

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def stitch(network, donor, P_network, P_donor, method='nearest',
           len_max=np.inf, label_suffix='', label_stitches='stitched'):
    r"""
    Stitches a second a network to the current network.

    Parameters
    ----------
    network : GenericNetwork
        The Network to which to donor Network will be attached
    donor : GenericNetwork
        The Network to stitch on to the current Network
    P_network : array_like
        The pores on the current Network
    P_donor : array_like
        The pores on the donor Network
    label_suffix : str or None
        Some text to append to each label in the donor Network before
        inserting them into the recipient.  The default is to append no
        text, but a common option would be to append the donor Network's
        name. To insert none of the donor labels, use ``None``.
    label_stitches : str or list[str]
        The label to apply to the newly created 'stitch' throats.  The
        defaul is 'stitched'.  If performing multiple stitches in a row it
        might be helpful to the throats created during each step uniquely
        for later identification.
    len_max : float
        Set a length limit on length of new throats
    method : str (default = 'nearest')
        The method to use when making pore to pore connections. Options are:

        - 'radius' : Connects each pore on the recipient network to the
                     nearest pores on the donor network, within ``len_max``
        - 'nearest' : Connects each pore on the recipienet network to the
                      nearest pore on the donor network.

    Notes
    -----
    Before stitching it is necessary to translate the pore coordinates of
    one of the Networks so that it is positioned correctly relative to the
    other.  This is illustrated in the example below.

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
    ...                     P_donor=pn2.pores('bottom'), method='radius',
    ...                     len_max=1.0)
    >>> [pn.Np, pn.Nt]
    [250, 625]

    """
    # Parse inputs
    if isinstance(label_stitches, str):
        label_stitches = [label_stitches]
    for s in label_stitches:
        if s not in network.keys():
            network['throat.' + s] = False
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
        [P1_ind, P2_ind] = np.where(D == D.min(axis=0))
        conns = np.vstack((P1[P1_ind], P2[P2_ind])).T
    elif method == 'radius':
        P1 = P_network
        P2 = P_donor + N_init['pore']  # Increment pores on donor
        C1 = network['pore.coords'][P_network]
        C2 = donor['pore.coords'][P_donor]
        D = sp.spatial.distance.cdist(C1, C2)
        [P1_ind, P2_ind] = np.where(D <= len_max)
        conns = np.vstack((P1[P1_ind], P2[P2_ind])).T
    else:
        raise Exception('<{}> method not supported'.format(method))

    merge_networks(network, donor)

    # Add the new stitch throats to the Network
    extend(network=network, throat_conns=conns, labels=label_stitches)


def stitch_pores(network, pores1, pores2, mode='gabriel'):
    r"""
    Stitches together pores in a network with disconnected clusters

    Parameters
    ----------
    network : OpenPNM Network
        The network to operate upon
    pores1 and pores2: array_like
        The pore indices of the disconnected clusters to be joined
    mode : str
        Dictates which tesselation method is used to identify which pores to
        stitch together.  Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'gabriel'    Uses the gabriel tesselation method
            'delaunay'   Uses the delaunay tesselation method
            ===========  =====================================================

    Returns
    -------
    None
        The network is operated on 'in-place' so nothing is returned.

    """
    from openpnm.network import Delaunay, Gabriel
    pores1 = network._parse_indices(pores1)
    pores2 = network._parse_indices(pores2)
    C1 = network.coords[pores1, :]
    C2 = network.coords[pores2, :]
    crds = np.vstack((C1, C2))
    if mode == 'delaunay':
        net = Delaunay(points=crds, settings={'trim': False})
    if mode == 'gabriel':
        net = Gabriel(points=crds, settings={'trim': False})
    net.set_label(pores=range(len(pores1)), label='pore.one')
    net.set_label(pores=range(len(pores2)), label='pore.two')
    Ts = net.find_neighbor_throats(pores=net.pores('one'), mode='xor')
    conns = net.conns[Ts]
    mapped_conns = np.vstack((pores1[conns[:, 0]],
                              pores2[conns[:, 1] - len(pores1)])).T
    mapped_conns = np.sort(mapped_conns, axis=1)
    extend(network=network, conns=mapped_conns, labels='stitched')


def connect_pores(network, pores1, pores2, labels=[], add_conns=True):
    r"""
    Returns the possible connections between two groups of pores, and optionally
    makes the connections.

    See ``Notes`` for advanced usage.

    Parameters
    ----------
    network : GenericNetwork

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
    (1) The method also works if ``pores1`` and ``pores2`` are list of lists,
    in which case it consecutively connects corresponding members of the two
    lists in a 1-to-1 fashion. Example: pores1 = [[0, 1], [2, 3]] and
    pores2 = [[5], [7, 9]] leads to creation of the following connections:

    ::

        0 --> 5     2 --> 7     3 --> 7
        1 --> 5     2 --> 9     3 --> 9

    (2) If you want to use the batch functionality, make sure that each element
    within ``pores1`` and ``pores2`` are of type list or ndarray.

    (3) It creates the connections in a format which is acceptable by
    the default OpenPNM connections ('throat.conns') and either adds them to
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

    """
    # Assert that `pores1` and `pores2` are list of lists
    try:
        len(pores1[0])
    except (TypeError, IndexError):
        pores1 = [pores1]
    try:
        len(pores2[0])
    except (TypeError, IndexError):
        pores2 = [pores2]

    if len(pores1) != len(pores2):
        raise Exception('Running in batch mode! pores1 and pores2 must be'
                        + ' of the same length.')

    arr1, arr2 = [], []
    for ps1, ps2 in zip(pores1, pores2):
        size1 = np.size(ps1)
        size2 = np.size(ps2)
        arr1.append(np.repeat(ps1, size2))
        arr2.append(np.tile(ps2, size1))
    conns = np.vstack([np.concatenate(arr1), np.concatenate(arr2)]).T
    if add_conns:
        extend(network=network, throat_conns=conns, labels=labels)
    else:
        return conns


def subdivide(network, pores, shape, labels=[]):
    r"""
    It trim the pores and replace them by cubic networks with the sent shape.

    Parameters
    ----------
    network : GenericNetwork
    pores : array_like
        The first group of pores to be replaced
    shape : array_like
        The shape of cubic networks in the target locations

    Notes
    -----
    It works only for cubic networks, and a check is performed to ensure this
    is the case.

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

    """
    mro = network._mro()
    if 'Cubic' not in mro:
        raise Exception('Subdivide is only supported for Cubic Networks')
    from openpnm.network import Cubic
    pores = network._parse_indices(pores)

    # Checks to find boundary pores in the selected pores
    if 'pore.boundary' in network.labels():
        if (np.in1d(pores, network.pores('boundary'))).any():
            raise Exception('boundary pores cannot be subdivided!')
    if not hasattr(network, '_subdivide_flag'):
        network._subdivide_flag = True
    else:
        raise Exception('The network has subdivided pores, so the method '
                        + 'does not support another subdivision')
    # Assigning right shape and division
    if np.size(shape) != 2 and np.size(shape) != 3:
        raise Exception('Subdivide not implemented for Networks other than 2D and 3D')
    if np.size(shape) == 3 and 1 not in shape:
        div = np.array(shape, ndmin=1)
        single_dim = None
    else:
        single_dim = np.where(np.array(get_shape(network)) == 1)[0]
        if np.size(single_dim) == 0:
            single_dim = None
        if np.size(shape) == 3:
            div = np.array(shape, ndmin=1)
        else:
            div = np.zeros(3, dtype=np.int32)
            if single_dim is None:
                dim = 2
            else:
                dim = single_dim
            div[dim] = 1
            div[-np.array(div, ndmin=1, dtype=bool)] = np.array(shape, ndmin=1)

    # Creating small network and handling labels
    networkspacing = get_spacing(network)
    new_netspacing = networkspacing/div
    new_net = Cubic(shape=div, spacing=new_netspacing)
    main_labels = ['front', 'back', 'left', 'right', 'top', 'bottom']
    if single_dim is not None:
        label_groups = np.array([['left', 'right'],
                                 ['front', 'back'],
                                 ['top', 'bottom']])
        non_single_labels = label_groups[np.array([0, 1, 2]) != single_dim]
    for label in main_labels:
        new_net['pore.surface_' + label] = False
        network['pore.surface_' + label] = False
        if single_dim is None:
            new_net['pore.surface_' + label][new_net.pores(labels=label)] = True
        else:
            for ind in [0, 1]:
                loc = (non_single_labels[ind] == label)
                temp_pores = new_net.pores(non_single_labels[ind][loc])
                new_net['pore.surface_' + label][temp_pores] = True

    old_coords = np.copy(new_net['pore.coords'])
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
        except KeyError:
            Pn_new_net = []
        Pn_old_net = Pn[~np.in1d(Pn, Pn_new_net)]
        Np1 = network.Np
        extend(pore_coords=new_net['pore.coords'],
               throat_conns=new_net['throat.conns'] + Np1,
               labels=labels, network=network)

        # Moving the temporary labels to the big network
        for label in main_labels:
            network['pore.surface_' + label][Np1:] = new_net['pore.surface_' + label]

        # Stitching the old pores of the main network to the new extended pores
        surf_pores = network.pores('surface_*')
        surf_coord = network['pore.coords'][surf_pores]
        for neighbor in Pn:
            neighbor_coord = network['pore.coords'][neighbor]
            dist = [round(np.inner(neighbor_coord-x, neighbor_coord-x),
                          20) for x in surf_coord]
            nearest_neighbor = surf_pores[dist == np.amin(dist)]
            if neighbor in Pn_old_net:
                coplanar_labels = network.labels(pores=nearest_neighbor)
                new_neighbors = network.pores(coplanar_labels,
                                              mode='and')
                # This might happen to the edge of the small network
                if np.size(new_neighbors) == 0:
                    labels = network.labels(pores=nearest_neighbor,
                                            mode='and')
                    common_label = [label for label in labels if 'surface_' in label]
                    new_neighbors = network.pores(common_label)
            elif neighbor in Pn_new_net:
                new_neighbors = nearest_neighbor
            connect_pores(network=network, pores1=neighbor,
                          pores2=new_neighbors, labels=labels)

        # Removing temporary labels
        for label in main_labels:
            network['pore.surface_' + label] = False
        new_net['pore.coords'] = np.copy(old_coords)

    label_faces(network=network)
    for label in main_labels:
        del network['pore.surface_' + label]
    trim(network=network, pores=pores)
    ws = network.project.workspace
    ws.close_project(new_net.project)


def merge_pores(network, pores, labels=['merged']):
    r"""
    Combines a selection of pores into a new single pore located at the
    centroid of the selected pores and connected to all of their neighbors.

    Parameters
    ----------
    network : GenericNetwork
    pores : array_like
        The list of pores which are to be combined into a new single pore
    labels : str or list[str]
        The labels to apply to the new pore and new throat connections

    Notes
    -----
    (1) The method also works if a list of lists is passed, in which case
    it consecutively merges the given selections of pores.

    (2) The selection of pores should be chosen carefully, preferrable so that
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
    # Assert that `pores` is list of lists
    try:
        len(pores[0])
    except (TypeError, IndexError):
        pores = [pores]

    N = len(pores)
    NBs, XYZs = [], []

    for Ps in pores:
        temp = network.find_neighbor_pores(pores=Ps, mode='union', flatten=True,
                                           include_input=False)
        NBs.append(temp)
        points = np.concatenate((temp, Ps))
        XYZs.append(hull_centroid(network["pore.coords"][points]))

    extend(network, pore_coords=XYZs, labels=labels)
    Pnew = network.Ps[-N::]

    # Possible throats between new pores: This only happens when running in
    # batch mode, i.e. multiple groups of pores are to be merged. In case
    # some of these groups share elements, possible throats between the
    # intersecting elements is not captured and must be added manually.
    pores_set = [set(items) for items in pores]
    NBs_set = [set(items) for items in NBs]
    ps1, ps2 = [], []
    from itertools import combinations
    for i, j in combinations(range(N), 2):
        if not NBs_set[i].isdisjoint(pores_set[j]):
            ps1.append([network.Ps[-N+i]])
            ps2.append([network.Ps[-N+j]])

    # Add (possible) connections between the new pores
    connect_pores(network, pores1=ps1, pores2=ps2, labels=labels)
    # Add connections between the new pores and the rest of the network
    connect_pores(network, pores2=np.split(Pnew, N), pores1=NBs, labels=labels)
    # Trim merged pores from the network
    trim(network=network, pores=np.concatenate(pores))


def add_reservoir_pore(cls, network, pores, offset=0.1):
    r"""
    Adds a single pore connected to all ``pores`` to act as a reservoir

    This function is mostly needed to make network compatible with the
    Statoil file format, which requires reservoir pores on the inlet and
    outlet faces.

    Parameters
    ----------
    network : GenericNetwork
        The network to which the reservoir pore should be added
    pores : array_like
        The pores to which the reservoir pore should be connected to
    offset : scalar
        Controls the distance which the reservoir is offset from the given
        ``pores``.  The total displacement is found from the network
        dimension normal to given ``pores``, multiplied by ``offset``.

    Returns
    -------
    project : list
        An OpenPNM project object with a new geometry object added to
        represent the newly added pore and throats.

    """
    from openpnm.geometry import GenericGeometry
    import openpnm.models.geometry as mods
    # Check if a label was given and fetch actual indices
    if isinstance(pores, str):
        # Convert 'face' into 'pore.face' if necessary
        if not pores.startswith('pore.'):
            pores = 'pore.' + pores
        pores = network.pores(pores)
    # Find coordinates of pores on given face
    coords = network['pore.coords'][pores]
    # Normalize the coordinates based on full network size
    c_norm = coords/network['pore.coords'].max(axis=0)
    # Identify axis of face by looking for dim with smallest delta
    diffs = np.amax(c_norm - np.average(c_norm, axis=0), axis=0)
    ax = np.where(diffs == diffs.min())[0][0]
    # Add new pore at center of domain
    new_coord = network['pore.coords'].mean(axis=0)
    domain_half_length = np.ptp(network['pore.coords'][:, ax])/2
    if coords[:, ax].mean() < network['pore.coords'][:, ax].mean():
        new_coord[ax] = new_coord[ax] - domain_half_length*(1 + offset)
    if coords[:, ax].mean() > network['pore.coords'][:, ax].mean():
        new_coord[ax] = new_coord[ax] + domain_half_length*(1 + offset)
    Ps = np.arange(network.Np, network.Np + 1)
    extend(network=network, coords=[new_coord], labels=['reservoir'])
    conns = [[P, network.Np-1] for P in pores]
    Ts = np.arange(network.Nt, network.Nt + len(conns))
    extend(network=network, conns=conns, labels=['reservoir'])
    # Compute the geometrical properties of the reservoir pore and throats
    # Confirm if network has any geometry props on it
    props = {'throat.length', 'pore.diameter', 'throat.volume'}
    if len(set(network.keys()).intersection(props)) > 0:
        raise Exception('Geometrical properties should be moved to a '
                        + 'geometry object first')
        # or just do this?:  geo = Imported(network=network)
    geo = GenericGeometry(network=network, pores=Ps, throats=Ts)
    geo.add_model(propname='pore.diameter',
                  model=mods.geometry.pore_size.largest_sphere)
    geo.add_model(propname='throat.diameter_temp',
                  model=mods.geometry.throat_size.from_neighbor_pores,
                  mode='min')
    geo.add_model(propname='throat.diameter',
                  model=mods.misc.scaled,
                  prop='throat.diameter_temp', factor=0.5)
    geo.add_model(propname='throat.volume',
                  model=mods.geometry.throat_volume.cylinder)
    return network.project



def add_boundary_pores(network, pores, offset=None, move_to=None,
                       apply_label='boundary'):
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
        be offset.  Either this, or ``move_to`` must be specified.
    move_to : 3 x 1 array
        The location to move the boundary pores to.  A value of ``None``
        indicates that no translation should be applied in that axis.  For
        instance, ``[None, None, 0]`` indicates that the boundary pores should
        moved along the z-axis to the specified location.  Either this or
        ``offset`` must be specified.
    apply_label : str
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
    Ps = np.array(pores, ndmin=1)
    if Ps.dtype is bool:
        Ps = network.to_indices(Ps)
    if np.size(pores) == 0:  # Handle an empty array if given
        return np.array([], dtype=np.int64)
    # Clone the specifed pores
    clone_pores(network=network, pores=Ps)
    newPs = network.pores('pore.clone')
    del network['pore.clone']
    newTs = network.throats('clone')
    del network['throat.clone']
    if offset is not None:  # Offset the cloned pores
        network['pore.coords'][newPs] += offset
    if move_to is not None:  # Move the cloned pores
        for i, d in enumerate(move_to):
            if d is not None:
                temp = network['pore.coords'][newPs]
                temp[:, i] = d
                network['pore.coords'][newPs] = temp
    # Apply labels to boundary pores (trim leading 'pores' if present)
    label = apply_label.split('.')[-1]
    plabel = 'pore.' + label
    tlabel = 'throat.' + label
    network[plabel] = False
    network[plabel][newPs] = True
    network[tlabel] = False
    network[tlabel][newTs] = True
