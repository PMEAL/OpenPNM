# -*- coding: utf-8 -*-
"""
===============================================================================
Network.tools.topology: Assorted topological manipulation methods
===============================================================================

"""
import scipy as _sp
import numpy as _np
import scipy.ndimage as _spim
from OpenPNM.Base import logging as _logging
from OpenPNM.Base import Workspace as _workspace
logger = _logging.getLogger(__name__)
_mgr = _workspace()


def extend(network, pore_coords=[], throat_conns=[], labels=[]):
    r'''
    Add individual pores and/or throats to the network from a list of coords
    or conns.  This is an in-place operation, meaning the received Network
    object will be altered directly.

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

    '''
    if (network._phases != []):
        raise Exception('Network has active Phases, cannot proceed')

    logger.info('Extending network')
    Np_old = network.num_pores()
    Nt_old = network.num_throats()
    Np = Np_old + int(_sp.size(pore_coords)/3)
    Nt = Nt_old + int(_sp.size(throat_conns)/2)
    # Adjust 'all' labels
    del network['pore.all'], network['throat.all']
    network['pore.all'] = _sp.ones((Np,), dtype=bool)
    network['throat.all'] = _sp.ones((Nt,), dtype=bool)
    # Add coords and conns
    if _sp.size(pore_coords) > 0:
        coords = _sp.vstack((network['pore.coords'], pore_coords))
        network['pore.coords'] = coords
    if _sp.size(throat_conns) > 0:
        conns = _sp.vstack((network['throat.conns'], throat_conns))
        network['throat.conns'] = conns
    # Increase size of any prop or label arrays on Network
    for item in list(network.keys()):
        if item.split('.')[1] not in ['coords', 'conns', 'all']:
            if item.split('.')[0] == 'pore':
                N = Np
            else:
                N = Nt
            if network[item].dtype == bool:
                temp = _sp.where(network[item])[0]
                network[item] = _sp.zeros((N,), dtype=bool)
                network[item][temp] = True
            elif network[item].dtype == object:
                temp = network[item]
                network[item] = _sp.ndarray((N,), dtype=object)
                network[item][_sp.arange(0, _sp.shape(temp)[0])] = temp
            else:
                temp = network[item]
                try:
                    network[item] = _sp.ones((N, _sp.shape(temp)[1]),
                                             dtype=float)*_sp.nan
                except:
                    network[item] = _sp.ones((N,), dtype=float)*_sp.nan
                network[item][_sp.arange(0, _sp.shape(temp)[0])] = temp
    # Apply labels, if supplied
    if labels != []:
        # Convert labels to list if necessary
        if type(labels) is str:
            labels = [labels]
        for label in labels:
            # Remove pore or throat from label, if present
            label = label.split('.')[-1]
            if _sp.size(pore_coords) > 0:
                Ps = _sp.r_[Np_old:Np]
                if 'pore.'+label not in network.labels():
                    network['pore.'+label] = False
                network['pore.'+label][Ps] = True
            if _sp.size(throat_conns) > 0:
                Ts = _sp.r_[Nt_old:Nt]
                if 'throat.'+label not in network.labels():
                    network['throat.'+label] = False
                network['throat.'+label][Ts] = True
    # Regnerate the adjacency matrices
    network._update_network()


def trim(network, pores=[], throats=[]):
    '''
    Remove pores or throats from the network.  This is an in-place operation,
    meaning the received Network object will be altered directly.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network from which pores or throats should be removed
    pores (or throats) : array_like
        A boolean mask of length Np (or Nt) or a list of indices of the
        pores (or throats) to be removed.

    Notes
    -----
    Trimming only adjusts Phase, Geometry, and Physics objects. Trimming a
    Network that has already been used to run simulations will break those
    simulation objects.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> pn.Np
    125
    >>> pn.Nt
    300
    >>> pn.trim(pores=[1])
    >>> pn.Np
    124
    >>> pn.Nt
    296

    '''
    mgr = network.workspace
    for net in mgr.networks():
        if net._parent is network:
            raise Exception('This Network has been cloned, cannot trim')
    if (_sp.size(pores) > 0) and (_sp.size(throats) > 0):
        raise Exception('Cannot delete pores and throats simultaneously')
    elif _sp.size(pores) > 0:
        pores = _sp.array(pores, ndmin=1)
        Pkeep = _sp.ones((network.num_pores(),), dtype=bool)
        Pkeep[pores] = False
        Tkeep = _sp.ones((network.num_throats(),), dtype=bool)
        Ts = network.find_neighbor_throats(pores)
        if len(Ts) > 0:
            Tkeep[Ts] = False
    elif _sp.size(throats) > 0:
        throats = _sp.array(throats, ndmin=1)
        Tkeep = _sp.ones((network.num_throats(),), dtype=bool)
        Tkeep[throats] = False
        Pkeep = network['pore.all'].copy()
    else:
        logger.warning('No pores or throats recieved')
        return

    # Trim all associated objects
    for item in network._geometries+network._physics+network._phases:
        Pnet = network['pore.'+item.name]*Pkeep
        Tnet = network['throat.'+item.name]*Tkeep
        temp = network.map_pores(pores=_sp.where(Pnet)[0],
                                 target=item,
                                 return_mapping=True)
        Ps = temp['target']
        temp = network.map_throats(throats=_sp.where(Tnet)[0],
                                   target=item,
                                   return_mapping=True)
        Ts = temp['target']
        # Then resize 'all
        item.update({'pore.all': _sp.ones((_sp.sum(Pnet),), dtype=bool)})
        item.update({'throat.all': _sp.ones((_sp.sum(Tnet),), dtype=bool)})
        # Overwrite remaining data and info
        for key in list(item.keys()):
            if key.split('.')[1] not in ['all']:
                temp = item.pop(key)
                if key.split('.')[0] == 'throat':
                    logger.debug('Trimming {a} from {b}'.format(a=key,
                                                                b=item.name))
                    item[key] = temp[Ts]
                if key.split('.')[0] == 'pore':
                    logger.debug('Trimming {a} from {b}'.format(a=key,
                                                                b=item.name))
                    item[key] = temp[Ps]

    # Remap throat connections
    Pmap = _sp.ones((network.Np,), dtype=int)*-1
    Pmap[Pkeep] = _sp.arange(0, _sp.sum(Pkeep))
    tpore1 = network['throat.conns'][:, 0]
    tpore2 = network['throat.conns'][:, 1]
    Tnew1 = Pmap[tpore1[Tkeep]]
    Tnew2 = Pmap[tpore2[Tkeep]]
    # Write 'all' label specifically
    network.update({'throat.all': _sp.ones((_sp.sum(Tkeep),), dtype=bool)})
    network.update({'pore.all': _sp.ones((_sp.sum(Pkeep),), dtype=bool)})
    # Write throat connections specifically
    network.update({'throat.conns': _sp.vstack((Tnew1, Tnew2)).T})
    # Overwrite remaining data and info
    for item in list(network.keys()):
        if item.split('.')[-1] not in ['conns', 'all']:
            temp = network.pop(item)
            if item.split('.')[0] == 'throat':
                logger.debug('Trimming {a} from {b}'.format(a=item,
                                                            b=network.name))
                network[item] = temp[Tkeep]
            if item.split('.')[0] == 'pore':
                logger.debug('Trimming {a} from {b}'.format(a=item,
                                                            b=network.name))
                network[item] = temp[Pkeep]

    # Reset network graphs
    network._update_network(mode='regenerate')

    # Check Network health
    health = network.check_network_health()
    if health['trim_pores'] != []:
        logger.warning('Isolated pores exist!  Run check_network_health to ID \
                        which pores to remove.')
        pass


def clone_pores(network, pores, apply_label=['clone'], mode='parents'):
    r'''
    Clones the specified pores and adds them to the network

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which the new pores are to be added
    pores : array_like
        List of pores to clone
    apply_labels : string, or list of strings
        The labels to apply to the clones, default is 'clone'
    mode : string
        Controls the connections between parents and clones.  Options are:

        - 'parents': (Default) Each clone is connected only to its parent
        - 'siblings': Clones are only connected to each other in the same
                      manner as parents were connected
        - 'isolated': No connections between parents or siblings
    '''
    if (network._geometries != []):
        logger.warning('Network has active Geometries, new pores must be \
                        assigned a Geometry')
    if (network._phases != []):
        raise Exception('Network has active Phases, cannot proceed')

    logger.debug('Cloning pores')
    apply_label = list(apply_label)
    # Clone pores
    Np = network.num_pores()
    Nt = network.num_throats()
    parents = _sp.array(pores, ndmin=1)
    pcurrent = network['pore.coords']
    pclone = pcurrent[pores, :]
    pnew = _sp.concatenate((pcurrent, pclone), axis=0)
    Npnew = _sp.shape(pnew)[0]
    clones = _sp.arange(Np, Npnew)
    # Add clone labels to network
    for item in apply_label:
        if 'pore.' + item not in network.keys():
            network['pore.'+item] = False
        if 'throat.' + item not in network.keys():
            network['throat.'+item] = False
    # Add connections between parents and clones
    if mode == 'parents':
        tclone = _sp.vstack((parents, clones)).T
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'siblings':
        ts = network.find_neighbor_throats(pores=pores, mode='intersection')
        tclone = network['throat.conns'][ts] + network.num_pores()
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'isolated':
        extend(network=network, pore_coords=pclone)
    # Apply provided labels to cloned pores
    for item in apply_label:
        network['pore.'+item][network.pores('all') >= Np] = True
        network['throat.'+item][network.throats('all') >= Nt] = True

    # Any existing adjacency and incidence matrices will be invalid
    network._update_network()


def stitch(network, donor, P_network, P_donor, method='nearest',
           len_max=_sp.inf, label_suffix=''):
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
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> pn2 = OpenPNM.Network.TestNet()
    >>> [pn.Np, pn.Nt]
    [125, 300]
    >>> [pn2.Np, pn2.Nt]
    [125, 300]
    >>> pn2['pore.coords'][:, 2] += 5.0
    >>> pn.stitch(donor=pn2, P_network=pn.pores('top'),
    ...           P_donor=pn2.pores('bottom'), method='nearest', len_max=1.0)
    >>> [pn.Np, pn.Nt]
    [250, 625]

    '''
    # Ensure Networks have no associated objects yet
    if (len(network._simulation()) > 1) or (len(donor._simulation()) > 1):
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
        D = _sp.spatial.distance.cdist(C1, C2)
        [P1_ind, P2_ind] = _sp.where(D <= len_max)
        conns = _sp.vstack((P1[P1_ind], P2[P2_ind])).T
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
    L = _sp.sum((C1 - C2)**2, axis=1)**0.5
    conns = conns[L <= len_max]

    # Add donor labels to recipient network
    if label_suffix is not None:
        if label_suffix != '':
            label_suffix = '_'+label_suffix
        for label in donor.labels():
            element = label.split('.')[0]
            locations = _sp.where(network._get_indices(element) >=
                                  N_init[element])[0]
            try:
                network[label + label_suffix]
            except:
                network[label + label_suffix] = False
            network[label+label_suffix][locations] = donor[label]

    # Add the new stitch throats to the Network
    extend(network=network, throat_conns=conns, labels='stitched')

    # Remove donor from Workspace, if present
    # This check allows for the reuse of a donor Network multiple times
    if donor in _mgr.values():
        _mgr.purge_object(donor)


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
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> pn.Nt
    300
    >>> pn.connect_pores(pores1=[22, 32], pores2=[16, 80, 68])
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
    size1 = _sp.size(pores1)
    size2 = _sp.size(pores2)
    array1 = _sp.repeat(pores1, size2)
    array2 = _sp.tile(pores2, size1)
    conns = _sp.vstack([array1, array2]).T
    if add_conns:
        extend(network=network, throat_conns=conns, labels=labels)
    else:
        return conns


def find_centroid(coords=None, mode='geometric'):
    r'''
    Finds the coordinates of the centroid of the specified pores.

    Parameters
    ----------
    coords : array_like
        An Np x 3 list of [x, y, z] coordinates for the pores of which the
        centroid is sought

    mode : string
        Controls how the centroid is calculated.  Options are:

        **'geometric'** : (default) Simply calcuates the average [X, Y, Z] of
        the given coordinates using scipy.mean(coords, axis=0)

        **'center_of_mass'** : not implemented yet

        **'inscribed_sphere'** : not implemented yet
    '''
    if mode == 'geometric':
        c = _sp.mean(coords, axis=0)
    return c


def find_pores_distance(network, pores1=None, pores2=None):
    r'''
    It finds the distance between two group of pores.
    '''
    from scipy.spatial.distance import cdist
    p1 = _sp.array(pores1, ndmin=1)
    p2 = _sp.array(pores2, ndmin=1)
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
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[5,6,5], spacing=0.001)
    >>> pn.Np
    150
    >>> nano_pores = [2,13,14,15]
    >>> pn.subdivide(pores=nano_pores, shape=[4,7,3], labels='nano')
    >>> pn.Np
    482
    >>> assert pn.Np == (150+4*(4*7*3)-4)

    '''
    mro = [item.__name__ for item in network.__class__.__mro__]
    if 'Cubic' not in mro:
        raise Exception('Subdivide is only supported for Cubic Networks')
    from OpenPNM.Network import Cubic
    pores = _sp.array(pores, ndmin=1)

    # Checks to find boundary pores in the selected pores
    if 'pore.boundary' in network.labels():
        if (_sp.in1d(pores, network.pores('boundary'))).any():
            raise Exception('boundary pores cannot be subdivided!')
    if not hasattr(network, '_subdivide_flag'):
        network._subdivide_flag = True
    else:
        raise Exception('The network has subdivided pores, so the method ' +
                        'does not support another subdivision.')
    # Assigning right shape and division
    if _sp.size(shape) != 2 and _sp.size(shape) != 3:
        raise Exception('Subdivide not implemented for Networks other than 2D \
                         and 3D')
    elif _sp.size(shape) == 3 and 1 not in shape:
        div = _sp.array(shape, ndmin=1)
        single_dim = None
    else:
        single_dim = _sp.where(_sp.array(network._shape) == 1)[0]
        if _sp.size(single_dim) == 0:
            single_dim = None
        if _sp.size(shape) == 3:
            div = _sp.array(shape, ndmin=1)
        else:
            div = _sp.zeros(3, dtype=_sp.int32)
            if single_dim is None:
                dim = 2
            else:
                dim = single_dim
            div[dim] = 1
            div[-_sp.array(div, ndmin=1, dtype=bool)] = _sp.array(shape,
                                                                  ndmin=1)

    # Creating small network and handling labels
    network_spacing = network._spacing
    new_net_spacing = network_spacing/div
    new_net = Cubic(shape=div, spacing=new_net_spacing)
    main_labels = ['left', 'right', 'front', 'back', 'top', 'bottom']
    if single_dim is not None:
        label_groups = _sp.array([['front', 'back'],
                                  ['left', 'right'],
                                  ['top', 'bottom']])
        non_single_labels = label_groups[_sp.array([0, 1, 2]) != single_dim]
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

    old_coords = _sp.copy(new_net['pore.coords'])
    if labels == []:
        labels = ['pore.subdivided_' + new_net.name]
    for P in pores:
        # Shifting the new network to the right location and attaching it to
        # the main network
        shift = network['pore.coords'][P] - network_spacing/2
        new_net['pore.coords'] += shift
        Pn = network.find_neighbor_pores(pores=P)
        try:
            Pn_new_net = network.pores(labels)
        except:
            Pn_new_net = []
        Pn_old_net = Pn[~_sp.in1d(Pn, Pn_new_net)]
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
            dist = [round(_sp.inner(neighbor_coord-x, neighbor_coord-x),
                          20) for x in surf_coord]
            nearest_neighbor = surf_pores[dist == _sp.amin(dist)]
            if neighbor in Pn_old_net:
                coplanar_labels = network.labels(pores=nearest_neighbor)
                new_neighbors = network.pores(coplanar_labels,
                                              mode='intersection')
                # This might happen to the edge of the small network
                if _sp.size(new_neighbors) == 0:
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
        new_net['pore.coords'] = _sp.copy(old_coords)

    network._label_surfaces()
    for l in main_labels:
        del network['pore.surface_'+l]
    trim(network=network, pores=pores)
    _mgr.purge_object(obj=new_net, mode='complete')


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
    if _sp.sum(occluded_ts) > 0:
        # Apply mask
        occluded_ts *= network["throat."+mask]
        trim(network=network, throats=occluded_ts)
    # Also get rid of isolated pores
    isolated_ps = network.check_network_health()['isolated_pores']
    if _sp.size(isolated_ps) > 0:
        # Convert to Bool array and apply mask
        temp_array = _sp.zeros(network.num_pores()).astype(bool)
        temp_array[isolated_ps] = True
        isolated_ps = temp_array * network["pore."+mask]
        trim(network=network, pores=isolated_ps)


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
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[20,20,1])
    >>> topo = op.Utilities.topology()
    >>> P = pn.find_nearby_pores(pores=111, distance=5, flatten=True)
    >>> topo.merge_pores(network=pn, pores=P, labels=['merged'])
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
    xyz = _sp.mean(network['pore.coords'][pores], axis=0)
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
    rmax = _np.array(outer_radius, ndmin=1)
    rmin = _np.array(inner_radius, ndmin=1)
    ind = 2 * rmax - 1
    coord = _np.indices((ind * _np.ones(dim, dtype=int)))
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
    return (img)


def template_sphere_shell(outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a sphere-shell. It is useful for
    passing to Cubic networks as a ``template`` to make spherical shaped
    networks.

    Parameters
    ----------
    outer_radius : int
        Number of the nodes in the outer radius of the shell

    inner_radius : int
        Number of the nodes in the inner radius of the shell

    Returns
    -------
    A Numpy array containing 1's to demarcate the sphere-shell, and 0's
    elsewhere.

    """
    img = _template_sphere_disc(dim=3, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return(img)


def template_disc_ring(outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a disc-ring.  It is useful for
    passing to Cubic networks as a ``template`` to make circular-shaped 2D
    networks.

    Parameters
    ----------
    outer_radius : int
        Number of the nodes in the outer radius of the disc

    inner_radius : int
        Number of the nodes in the inner radius of the disc

    Returns
    -------
    A Numpy array containing 1's to demarcate the disc-ring, and 0's
    elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return(img)


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
    >>> import OpenPNM as op
    >>> net = op.Network.Cubic(shape=[5, 5, 5])
    >>> op.Network.tools.find_surface_pores(network=net)
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
        (xmax, ymax, zmax) = _sp.amax(network['pore.coords'], axis=0)
        (xmin, ymin, zmin) = _sp.amin(network['pore.coords'], axis=0)
        xave = (xmin+xmax)/2
        yave = (ymin+ymax)/2
        zave = (zmin+zmax)/2
        markers = [[xmax + xave, yave, zave],
                   [xmin - xave, yave, zave],
                   [xave, ymax + yave, zave],
                   [xave, ymin - yave, zave],
                   [xave, yave, zmax + zave],
                   [xave, yave, zmin - zave]]
    markers = _sp.atleast_2d(markers)
    tri = sptl.Delaunay(network['pore.coords'], incremental=True)
    tri.add_points(markers)
    (indices, indptr) = tri.vertex_neighbor_vertices
    for k in range(network.Np, tri.npoints):
        neighbors = indptr[indices[k]:indices[k+1]]
        inds = _sp.where(neighbors < network.Np)
        neighbors = neighbors[inds]
        if 'pore.'+label not in network.keys():
            network['pore.'+label] = False
        network['pore.'+label][neighbors] = True


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
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundaries()
    >>> Ts = pn.throats('*boundary', mode='not')
    >>> # Create figure showing boundary throats
    >>> fig = op.Network.tools.plot_connections(network=pn, throats=Ts)
    >>> Ts = pn.throats('*boundary')
    >>> # Pass existing fig back into function to plot additional throats
    >>> fig = op.Network.tools.plot_connections(network=pn, throats=Ts,
    ...                                         fig=fig, color='r')

    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if throats is None:
        Ts = network.Ts
    else:
        Ts = network._parse_locations(locations=throats)

    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.get_axes()[0]

    # Create dummy indexing to sp.inf
    i = -1*_sp.ones((_sp.size(Ts)*3, ), dtype=int)
    i[0::3] = network['throat.conns'][Ts, 0]
    i[1::3] = network['throat.conns'][Ts, 1]

    # Collect coordinates and scale axes to fit
    Ps = _sp.unique(network['throat.conns'][Ts])
    X = network['pore.coords'][Ps, 0]
    Y = network['pore.coords'][Ps, 1]
    Z = network['pore.coords'][Ps, 2]
    _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)

    # Add sp.inf to the last element of pore.coords (i.e. -1)
    inf = _sp.array((_sp.inf,))
    X = _sp.hstack([network['pore.coords'][:, 0], inf])
    Y = _sp.hstack([network['pore.coords'][:, 1], inf])
    Z = _sp.hstack([network['pore.coords'][:, 2], inf])
    ax.plot(xs=X[i], ys=Y[i], zs=Z[i], **kwargs)

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
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 3])
    >>> pn.add_boundaries()
    >>> Ps = pn.pores('internal')
    >>> # Create figure showing internal pores
    >>> fig = op.Network.tools.plot_coordinates(network=pn, pores=Ps, c='b')
    >>> Ps = pn.pores('*boundary')
    >>> # Pass existing fig back into function to plot boundary pores
    >>> fig = op.Network.tools.plot_coordinates(network=pn, pores=Ps, fig=fig,
    ...                                         c='r')

    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if pores is None:
        Ps = network.Ps
    else:
        Ps = network._parse_locations(locations=pores)

    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.get_axes()[0]

    # Collect specified coordinates
    X = network['pore.coords'][Ps, 0]
    Y = network['pore.coords'][Ps, 1]
    Z = network['pore.coords'][Ps, 2]
    _scale_3d_axes(ax=ax, X=X, Y=Y, Z=Z)

    ax.scatter(xs=X, ys=Y, zs=Z, **kwargs)
    return fig


def _scale_3d_axes(ax, X, Y, Z):
    if hasattr(ax, '_scaled'):
        logger.warning('Axes is already scaled to previously plotted data')
    else:
        ax._scaled = True
        max_range = _sp.array([X.max()-X.min(), Y.max()-Y.min(),
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

    >>> import OpenPNM as op
    >>> import scipy as sp
    >>> import scipy.ndimage as spim
    >>> im = sp.ones([21, 21, 21], dtype=int)
    >>> im[10, 10, 10] = 0
    >>> im = spim.distance_transform_edt(im) <= 20  # Create sphere of 1's
    >>> prob = spim.distance_transform_edt(im)
    >>> prob = prob / sp.amax(prob)  # Normalize between 0 and 1
    >>> pts = op.Network.tools.generate_base_points(num_points=50,
    ...                                             domain_size=[2],
    ...                                             prob=prob)
    >>> net = op.Network.DelaunayVoronoiDual(points=pts, domain_size=[2])
    """
    def _try_points(num_points, prob):
        prob = _sp.array(prob)/_sp.amax(prob)  # Ensure prob is normalized
        base_pts = []
        N = 0
        while N < num_points:
            pt = _sp.random.rand(3)  # Generate a point
            # Test whether to keep it or not
            [indx, indy, indz] = _sp.floor(pt*_sp.shape(prob)).astype(int)
            if _sp.random.rand(1) <= prob[indx][indy][indz]:
                base_pts.append(pt)
                N += 1
        base_pts = _sp.array(base_pts)
        return base_pts
    if len(domain_size) == 1:  # Spherical
        domain_size = _sp.array(domain_size)
        if prob is None:
            prob = _sp.ones([41, 41, 41])
            prob[20, 20, 20] = 0
            prob = _spim.distance_transform_bf(prob) <= 20
        base_pts = _try_points(num_points, prob)
        # Convert to spherical coordinates
        [X, Y, Z] = _sp.array(base_pts - [0.5, 0.5, 0.5]).T  # Center at origin
        r = 2*_sp.sqrt(X**2 + Y**2 + Z**2)*domain_size[0]
        theta = 2*_sp.arctan(Y/X)
        phi = 2*_sp.arctan(_sp.sqrt(X**2 + Y**2)/Z)
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, phi] = [r[inds], theta[inds], phi[inds]]
        # Reflect base points across perimeter
        new_r = 2*domain_size - r
        r = _sp.hstack([r, new_r])
        theta = _sp.hstack([theta, theta])
        phi = _sp.hstack([phi, phi])
        # Convert to Cartesean coordinates
        X = r*_sp.cos(theta)*_sp.sin(phi)
        Y = r*_sp.sin(theta)*_sp.sin(phi)
        Z = r*_sp.cos(phi)
        base_pts = _sp.vstack([X, Y, Z]).T
    elif len(domain_size) == 2:  # Cylindrical
        domain_size = _sp.array(domain_size)
        if prob is None:
            prob = _sp.ones([41, 41, 41])
            prob[20, 20, :] = 0
            prob = _spim.distance_transform_bf(prob) <= 20
        base_pts = _try_points(num_points, prob)
        # Convert to cylindrical coordinates
        [X, Y, Z] = _sp.array(base_pts - [0.5, 0.5, 0]).T  # Center on z-axis
        r = 2*_sp.sqrt(X**2 + Y**2)*domain_size[0]
        theta = 2*_sp.arctan(Y/X)
        z = Z*domain_size[1]
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        inds = ~((z > domain_size[1]) + (z < 0))
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        # Reflect base points about faces and perimeter
        new_r = 2*domain_size[0] - r
        r = _sp.hstack([r, new_r])
        theta = _sp.hstack([theta, theta])
        z = _sp.hstack([z, z])
        r = _sp.hstack([r, r, r])
        theta = _sp.hstack([theta, theta, theta])
        z = _sp.hstack([z, -z, 2-z])
        # Convert to Cartesean coordinates
        X = r*_sp.cos(theta)
        Y = r*_sp.sin(theta)
        Z = z
        base_pts = _sp.vstack([X, Y, Z]).T
    elif len(domain_size) == 3:  # Rectilinear
        domain_size = _sp.array(domain_size)
        Nx, Ny, Nz = domain_size
        if prob is None:
            prob = _sp.ones([10, 10, 10], dtype=float)
        base_pts = _try_points(num_points, prob)
        base_pts = base_pts*domain_size
        # Reflect base points about all 6 faces
        orig_pts = base_pts
        base_pts = _sp.vstack((base_pts, [-1, 1, 1]*orig_pts +
                                         [2.0*Nx, 0, 0]))
        base_pts = _sp.vstack((base_pts, [1, -1, 1]*orig_pts +
                                         [0, 2.0*Ny, 0]))
        base_pts = _sp.vstack((base_pts, [1, 1, -1]*orig_pts +
                                         [0, 0, 2.0*Nz]))
        base_pts = _sp.vstack((base_pts, [-1, 1, 1]*orig_pts))
        base_pts = _sp.vstack((base_pts, [1, -1, 1]*orig_pts))
        base_pts = _sp.vstack((base_pts, [1, 1, -1]*orig_pts))
    return base_pts
