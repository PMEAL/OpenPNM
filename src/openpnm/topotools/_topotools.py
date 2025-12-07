import logging
import numpy as np
import scipy as sp
from scipy.spatial import cKDTree
from scipy.sparse import csgraph
from scipy.spatial import ConvexHull
from openpnm.utils import Workspace
import openpnm._skgraph as skgr


logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = [
    'isoutside',
    'rotate_coords',
    'shear_coords',
    'trim',
    'extend',
    'label_faces',
    'find_surface_pores',
    'dimensionality',
    'clone_pores',
    'merge_networks',
    'stitch',
    'connect_pores',
    'merge_pores',
    'hull_centroid',
    'template_sphere_shell',
    'template_cylinder_annulus',
    'generate_base_points',
    'reflect_base_points',
    'add_boundary_pores',
    'iscoplanar',
    'is_fully_connected',
    'get_spacing',
    'get_shape',
    'get_domain_area',
    'get_domain_length',
    'filter_pores_by_z',
    'find_interface_throats',
    'add_reservoir_pore',
    'reduce_coordination',
]


def isoutside(**kwargs):
    return skgr.tools.isoutside(**kwargs)


isoutside.__doc__ = skgr.tools.isoutside.__doc__


def rotate_coords(network, **kwargs):
    network['pore.coords'] = skgr.tools.rotate_coords(network.coords, **kwargs)
    return network


rotate_coords.__doc__ = skgr.tools.rotate_coords.__doc__


def shear_coords(network, **kwargs):
    network['pore.coords'] = skgr.tools.shear_coords(network.coords, **kwargs)
    return network


shear_coords.__doc__ = skgr.tools.shear_coords.__doc__


def template_sphere_shell(r_outer, r_inner=0):
    return skgr.generators.tools.template_sphere_shell(r_outer, r_inner)


template_sphere_shell.__doc__ = \
    skgr.generators.tools.template_sphere_shell.__doc__


def template_cylinder_annulus(z, r_outer, r_inner=0):
    return skgr.generators.tools.template_cylinder_annulus(z, r_outer, r_inner)


template_cylinder_annulus.__doc__ = \
    skgr.generators.tools.template_cylinder_annulus.__doc__


def generate_base_points(num_points, domain_size, reflect=True):
    return skgr.generators.tools.generate_base_points(
        num_points,
        domain_size,
        reflect,
    )


generate_base_points.__doc__ = \
    skgr.generators.tools.generate_base_points.__doc__


def reflect_base_points(points, domain_size):
    return skgr.generators.tools.reflect_base_points(points, domain_size)


reflect_base_points.__doc__ = \
    skgr.generators.tools.reflect_base_points.__doc__


def get_spacing(network):
    return skgr.tools.get_cubic_spacing(network)


get_spacing.__doc__ = skgr.tools.get_cubic_spacing.__doc__


def get_shape(network):
    return skgr.tools.get_cubic_shape(network)


get_shape.__doc__ = skgr.tools.get_cubic_shape.__doc__


def filter_pores_by_z(network, pores, z):
    return skgr.queries.filter_by_z(network=network, inds=pores, z=z)


filter_pores_by_z.__doc__ = skgr.queries.filter_by_z.__doc__


def find_interface_throats(network, P1, P2):
    return skgr.queries.find_common_edges(network=network, inds_1=P1, inds_2=P2)


find_interface_throats.__doc__ = skgr.queries.find_common_edges.__doc__


def dimensionality(network):
    r"""
    Determines whether a network is 1D, 2D or 3D, and in which dimensions

    Parameters
    ----------
    network : dict
        The OpenPNM network of interest

    Returns
    -------
    dims : boolean mask
        A 3 x 1 array of booleans with ``True`` indicating if any
        dimensionality exists on each axis.

    Notes
    -----
    This function looks at the coordinates of each pore and if any axes all
    have the same values that axis is considered non-dimensional.

    """
    return skgr.tools.dimensionality(network)


def trim(network, pores=[], throats=[]):
    """
    Remove pores or throats from the network

    Parameters
    ----------
    network : Network
        The Network from which pores or throats should be removed
    pores (or throats) : array_like
        The indices of the of the pores or throats to be removed from the
        network.

    """
    pores = network._parse_indices(pores)
    throats = network._parse_indices(throats)
    Pkeep = np.copy(network['pore.all'])
    Tkeep = np.copy(network['throat.all'])
    if np.size(pores) > 0:
        Pkeep[pores] = False
        if not np.any(Pkeep):
            raise Exception('Cannot delete ALL pores')
        # # Performing customized find_neighbor_throats which is much faster,
        # # but not general for other types of queries
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
            logger.info('Removing ALL throats from network')
            for item in list(network.keys()):
                if item.split('.', 1)[0] == 'throat':
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
    for obj in network.project:
        if (obj.Np == Np_old) and (obj.Nt == Nt_old):
            Ps = Pkeep_inds
            Ts = Tkeep_inds
        for key in list(obj.keys()):
            temp = obj.pop(key)
            if key.split('.', 1)[0] == 'throat':
                obj.update({key: temp[Ts]})
            if key.split('.', 1)[0] == 'pore':
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
    network : Network
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
    objs = list(network.project.phases)
    objs.append(network)
    for obj in objs:
        obj.update({'pore.all': np.ones([Np, ], dtype=bool),
                    'throat.all': np.ones([Nt, ], dtype=bool)})
        for item in list(obj.keys()):
            N = obj._count(element=item.split('.', 1)[0])
            if obj[item].shape[0] < N:
                arr = obj.pop(item)
                s = arr.shape
                if arr.dtype == bool:
                    obj[item] = np.zeros(shape=(N, *s[1:]), dtype=bool)
                else:
                    obj[item] = np.ones(shape=(N, *s[1:]), dtype=float)*np.nan
                obj[item][:arr.shape[0]] = arr

    # Regenerate models on all objects to fill new elements
    for obj in network.project.phases:
        if hasattr(obj, 'models'):
            obj.regenerate_models()

    # Apply labels, if supplied
    if labels != []:
        # Convert labels to list if necessary
        if isinstance(labels, str):
            labels = [labels]
        for label in labels:
            # Remove pore or throat from label, if present
            label = label.split('.', 1)[-1]
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


def label_faces(network, tol=0.0, label='surface'):
    r"""
    Finds pores on the surface of the network and labels them according to
    whether they are on the *top*, *bottom*, etc.

    This function assumes the network is cubic in shape

    Parameters
    ----------
    network : Network
        The network to apply the labels
    tol : scalar
        The tolerance for defining what counts as a surface pore, which is
        specifically meant for random networks.  All pores with ``tol`` of
        the maximum or minimum along each axis are counts as pores.  The
        default is 0.
    label : str
        An identifying label to isolate the pores on the faces of the network.
        The default is 'surface'.  Surface pores can be found using
        ``find_surface_pores``.

    """
    if label is not None:
        label = label.split('.', 1)[-1]
        if 'pore.'+label not in network.labels():
            find_surface_pores(network, label=label)
        Psurf = network['pore.'+label]
    else:
        Psurf = True  # So it will "do nothing" below
    crds = network['pore.coords']
    xmin, xmax = np.amin(crds[:, 0]), np.amax(crds[:, 0])
    xspan = xmax - xmin
    ymin, ymax = np.amin(crds[:, 1]), np.amax(crds[:, 1])
    yspan = ymax - ymin
    zmin, zmax = np.amin(crds[:, 2]), np.amax(crds[:, 2])
    zspan = zmax - zmin
    dims = dimensionality(network)
    if dims[0]:
        network['pore.left'] = (crds[:, 0] <= (xmin + tol*xspan)) * Psurf
        network['pore.xmin'] = np.copy(network['pore.left'])
        network['pore.right'] = (crds[:, 0] >= (xmax - tol*xspan)) * Psurf
        network['pore.xmax'] = np.copy(network['pore.right'])
    if dims[1]:
        network['pore.front'] = (crds[:, 1] <= (ymin + tol*yspan)) * Psurf
        network['pore.ymin'] = np.copy(network['pore.front'])
        network['pore.back'] = (crds[:, 1] >= (ymax - tol*yspan)) * Psurf
        network['pore.ymax'] = np.copy(network['pore.back'])
    if dims[2]:
        network['pore.bottom'] = (crds[:, 2] <= (zmin + tol*zspan)) * Psurf
        network['pore.zmin'] = np.copy(network['pore.bottom'])
        network['pore.top'] = (crds[:, 2] >= (zmax - tol*zspan)) * Psurf
        network['pore.zmax'] = np.copy(network['pore.top'])


def find_surface_pores(network, markers=None, label='surface'):
    r"""
    Find the pores on the surface of the domain by performing a Delaunay
    triangulation between the network pores and some external ``markers``. All
    pores connected to these external marker points are considered surface
    pores.

    Parameters
    ----------
    network: Network
        The network for which the surface pores are to be found
    markers: array_like
        3 x N array of the marker coordinates to use in the triangulation.  The
        labeling is performed in one step, so all points are added, and then
        any pores connected to at least one marker is given the provided label.
        By default, this function will automatically generate 6 points outside
        each axis of the network domain. Users may wish to specify a single
        external marker point and provide an appropriate label in order to
        identify specific faces.  For instance, the marker may be *above* the
        domain, and the label might be 'top_surface'.
    label : str
        The label to apply to the pores.  The default is 'surface'.

    Notes
    -----
    This function does not check whether the given markers actually lie outside
    the domain, allowing the labeling of *internal* sufaces.

    If this method fails to mark some surface pores, consider sending more
    markers on each face.

    """
    import scipy.spatial as sptl
    dims = dimensionality(network)
    coords = network['pore.coords'][:, dims]
    if markers is None:
        # normalize coords to a 1 unit cube centered on origin
        coords -= np.amin(coords, axis=0)
        coords /= np.amax(coords, axis=0)
        coords -= 0.5
        npts = max((network.Np/10, 100))
        if sum(dims) == 1:
            network['pore.'+label] = True
            return
        if sum(dims) == 2:
            r = 0.75
            theta = np.linspace(0, 2*np.pi, int(npts), dtype=float)
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            markers = np.vstack((x, y)).T
        if sum(dims) == 3:
            r = 1.00
            indices = np.arange(0, int(npts), dtype=float) + 0.5
            phi = np.arccos(1 - 2*indices/npts)
            theta = np.pi * (1 + 5**0.5) * indices
            x = r*np.cos(theta) * np.sin(phi)
            y = r*np.sin(theta) * np.sin(phi)
            z = r*np.cos(phi)
            markers = np.vstack((x, y, z)).T
    else:
        if sum(dims) == 1:
            pass
        if sum(dims) == 2:
            markers = np.atleast_2d(markers)
            if markers.shape[1] != 2:
                raise Exception('Network appears planar, so markers must be 2D')
        if sum(dims) == 3:
            markers = np.atleast_2d(markers)
            if markers.shape[1] != 3:
                raise Exception('Markers must be 3D for this network')
    pts = np.vstack((coords, markers))
    tri = sptl.Delaunay(pts, incremental=False)
    (indices, indptr) = tri.vertex_neighbor_vertices
    for k in range(network.Np, tri.npoints):
        neighbors = indptr[indices[k]:indices[k+1]]
        inds = np.where(neighbors < network.Np)
        neighbors = neighbors[inds]
        if 'pore.'+label not in network.keys():
            network['pore.'+label] = False
        network['pore.'+label][neighbors] = True


def clone_pores(network, pores, labels=['clone'], mode='parents'):
    r"""
    Clones the specified pores and adds them to the network

    Parameters
    ----------
    network : Network
        The Network object to which the new pores are to be added
    pores : array_like
        List of pores to clone
    labels : str, or list[str]
        The labels to apply to the clones, default is 'clone'
    mode : str
        Controls the connections between parents and clones.  Options are:

        ===========  ==========================================================
        mode         description
        ===========  ==========================================================
        'parents'    Each clone is connected only to its parent.(Default)
        'siblings'   Clones are only connected to each other in the same
                     manner as parents were connected.
        'isolated'   No connections between parents or siblings
        ===========  ==========================================================

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
    Combine multiple networks into one without making any topological
    connections

    Parameters
    ----------
    network : Network
        The network to which all the other networks should be added.
    donor : Network or list of Objects
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

    for donor in donors:
        network['throat.conns'] = np.vstack((network['throat.conns'],
                                             donor['throat.conns']
                                             + network.Np))
        network['pore.coords'] = np.vstack((network['pore.coords'],
                                            donor['pore.coords']))
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
                        logger.debug('Adding ' + key + ' to donor')
                        if network[key].dtype == bool:  # Deal with labels
                            donor[key] = False
                        else:  # Deal with numerical data
                            element = key.split('.', 1)[0]
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
                    logger.debug('Adding ' + key + ' to network')
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
    network : Network
        The Network to which to donor Network will be attached
    donor : Network
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

        ========= =============================================================
        mode      description
        ========= =============================================================
        'radius'  Connects each pore on the recipient network to the nearest
                  pores on the donor network, within ``len_max``
        'nearest' Connects each pore on the recipienet network to the nearest
                  pore on the donor network.
        ========= =============================================================

    Notes
    -----
    Before stitching it is necessary to translate the pore coordinates of
    one of the Networks so that it is positioned correctly relative to the
    other.  This is illustrated in the example below.

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

    # Remove donor from Workspace, if present
    # This check allows for the reuse of a donor Network multiple times
    for sim in list(ws.values()):
        if donor in sim:
            del ws[sim.name]


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

        ===========  ==========================================================
        mode         meaning
        ===========  ==========================================================
        'delaunay'   Uses the delaunay tesselation method
        ===========  ==========================================================

    Returns
    -------
    None
        The network is operated on 'in-place' so nothing is returned.

    """
    raise NotImplementedError()
    from openpnm.network import Delaunay
    pores1 = network._parse_indices(pores1)
    pores2 = network._parse_indices(pores2)
    C1 = network.coords[pores1, :]
    C2 = network.coords[pores2, :]
    crds = np.vstack((C1, C2))
    if mode == 'delaunay':
        shape = get_shape(network) * dimensionality(network)
        net = Delaunay(points=crds, shape=shape)
    else:
        raise Exception('Unsupported mode')
    net.set_label(pores=range(len(pores1)), label='pore.one')
    net.set_label(pores=range(len(pores2)), label='pore.two')
    Ts = net.find_neighbor_throats(pores=net.pores('one'), mode='xor')
    conns = net.conns[Ts]
    mapped_conns = np.vstack((pores1[conns[:, 0]],
                              pores2[conns[:, 1] - len(pores1)])).T
    mapped_conns = np.sort(mapped_conns, axis=1)
    extend(network=network, conns=mapped_conns, labels='stitched')


def connect_pores(network, pores1, pores2, labels=['new_conns']):
    r"""
    Returns the possible connections between two groups of pores

    Parameters
    ----------
    network : Network
        The network to which the pores should be added
    pores1 : array_like
        The first group of pores on the network
    pores2 : array_like
        The second group of pores on the network
    labels : list of strings
        The labels to apply to the new throats. The default is ``'new_conns'``.

    Notes
    -----
    The method also works if ``pores1`` and ``pores2`` are list of lists,
    in which case it consecutively connects corresponding members of the two
    lists in a 1-to-1 fashion. Example: pores1 = [[0, 1], [2, 3]] and
    pores2 = [[5], [7, 9]] leads to creation of the following connections:

    ::

        0 --> 5     2 --> 7     3 --> 7
        1 --> 5     2 --> 9     3 --> 9

    If you want to use the batch functionality, make sure that each element
    within ``pores1`` and ``pores2`` are of type list or ndarray.

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
    extend(network=network, throat_conns=conns, labels=labels)


def merge_pores(network, pores, labels=['merged'], include_neighbors=True):
    r"""
    Combines a selection of pores into a new single pore located at the
    centroid of the selected pores (and optionally their neighbors)
    and connected to all of their neighbors.

    Parameters
    ----------
    network : Network
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

    """
    # Assert that `pores` is list of lists
    try:
        len(pores[0])
    except (TypeError, IndexError):
        pores = [pores]

    N = len(pores)
    NBs, XYZs = [], []

    for Ps in pores:
        temp = network.find_neighbor_pores(pores=Ps,
                                           mode='union',
                                           flatten=True,
                                           include_input=False)
        NBs.append(temp)
        if len(Ps) == 2:
            XYZs.append(np.mean(network["pore.coords"][Ps], axis=0))
        else:
            if include_neighbors:
                points = np.concatenate((temp, Ps))
            else:
                points = Ps
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


def hull_centroid(points):
    r"""
    Computes centroid of the convex hull enclosing the given coordinates.

    Parameters
    ----------
    points : Np by 3 ndarray
        Coordinates (xyz)

    Returns
    -------
    centroid : array
        A 3 by 1 Numpy array containing coordinates of the centroid.

    """
    dim = [np.unique(points[:, i]).size != 1 for i in range(3)]
    hull = ConvexHull(points[:, dim])
    centroid = points.mean(axis=0)
    centroid[dim] = hull.points[hull.vertices].mean(axis=0)

    return centroid


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
    label = apply_label.split('.', 1)[-1]
    plabel = 'pore.' + label
    tlabel = 'throat.' + label
    network[plabel] = False
    network[plabel][newPs] = True
    network[tlabel] = False
    network[tlabel][newTs] = True


def iscoplanar(coords):
    r"""
    Determines if given pores are coplanar with each other

    Parameters
    ----------
    coords : array_like
        List of pore coords to check for coplanarity.  At least 3 pores are
        required.

    Returns
    -------
    results : bool
        A boolean value of whether given points are coplanar (``True``) or
        not (``False``)

    """
    coords = np.array(coords, ndmin=1)
    if np.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:, 0]
    Py = coords[:, 1]
    Pz = coords[:, 2]

    # Do easy check first, for common coordinate
    if np.shape(np.unique(Px))[0] == 1:
        return True
    if np.shape(np.unique(Py))[0] == 1:
        return True
    if np.shape(np.unique(Pz))[0] == 1:
        return True

    # Perform rigorous check using vector algebra
    # Grab first basis vector from list of coords
    n1 = np.array((Px[1] - Px[0], Py[1] - Py[0], Pz[1] - Pz[0])).T
    n = np.array([0.0, 0.0, 0.0])
    i = 1
    while n.sum() == 0:
        if i >= (np.size(Px) - 1):
            logger.warning('No valid basis vectors found')
            return False
        # Chose a secon basis vector
        n2 = np.array((Px[i+1] - Px[i], Py[i+1] - Py[i], Pz[i+1] - Pz[i])).T
        # Find their cross product
        n = np.cross(n1, n2)
        i += 1
    # Create vectors between all other pairs of points
    r = np.array((Px[1:-1] - Px[0], Py[1:-1] - Py[0], Pz[1:-1] - Pz[0]))
    # Ensure they all lie on the same plane
    n_dot = np.dot(n, r)

    return bool(np.sum(np.absolute(n_dot)) == 0)


def is_fully_connected(network, pores_BC=None):
    r"""
    Checks whether network is fully connected, i.e. not clustered.

    Parameters
    ----------
    network : Network
        The network whose connectivity to check.
    pores_BC : array_like (optional)
        The pore indices of boundary conditions (inlets/outlets).

    Returns
    -------
    bool
        If ``pores_BC`` is not specified, then returns ``True`` only if
        the entire network is connected to the same cluster. If
        ``pores_BC`` is given, then returns ``True`` only if all clusters
        are connected to the given boundary condition pores.

    """
    am = network.get_adjacency_matrix(fmt='lil').copy()
    temp = csgraph.connected_components(am, directed=False)[1]
    is_connected = np.unique(temp).size == 1
    # Ensure all clusters are part of pores, if given
    if not is_connected and pores_BC is not None:
        am.resize(network.Np + 1, network.Np + 1)
        pores_BC = network._parse_indices(pores_BC)
        am.rows[-1] = pores_BC.tolist()
        am.data[-1] = np.arange(network.Nt, network.Nt + len(pores_BC)).tolist()
        temp = csgraph.connected_components(am, directed=False)[1]
        is_connected = np.unique(temp).size == 1
    return is_connected


def get_domain_area(network, inlets=None, outlets=None):
    r"""
    Determine the cross sectional area relative to the inlets/outlets.

    Parameters
    ----------
    network : Network
        The network object containing the pore coordinates
    inlets : array_like
        The pore indices of the inlets.
    outlets : array_Like
        The pore indices of the outlets.

    Returns
    -------
    area : scalar
        The cross sectional area relative to the inlets/outlets.

    """
    logger.warning('Attempting to estimate inlet area...will be low')
    if dimensionality(network).sum() != 3:
        raise Exception('The network is not 3D, specify area manually')
    inlets = network.coords[inlets]
    outlets = network.coords[outlets]
    if not iscoplanar(inlets):
        logger.error('Detected inlet pores are not coplanar')
    if not iscoplanar(outlets):
        logger.error('Detected outlet pores are not coplanar')
    Nin = np.ptp(inlets, axis=0) > 0
    if Nin.all():
        logger.warning('Detected inlets are not oriented along a principle axis')
    Nout = np.ptp(outlets, axis=0) > 0
    if Nout.all():
        logger.warning('Detected outlets are not oriented along a principle axis')
    hull_in = ConvexHull(points=inlets[:, Nin])
    hull_out = ConvexHull(points=outlets[:, Nout])
    if hull_in.volume != hull_out.volume:
        logger.error('Inlet and outlet faces are different area')
    area = hull_in.volume  # In 2D: volume=area, area=perimeter
    return area


def get_domain_length(network, inlets=None, outlets=None):
    r"""
    Determine the domain length relative to the inlets/outlets.

    Parameters
    ----------
    network : Network
        The network object containing the pore coordinates
    inlets : array_like
        The pore indices of the inlets.
    outlets : array_Like
        The pore indices of the outlets.

    Returns
    -------
    area : scalar
        The domain length relative to the inlets/outlets.

    """
    msg = ('Attempting to estimate domain length...could be low if'
           ' boundary pores were not added')
    logger.warning(msg)
    inlets = network.coords[inlets]
    outlets = network.coords[outlets]
    if not iscoplanar(inlets):
        logger.error('Detected inlet pores are not coplanar')
    if not iscoplanar(outlets):
        logger.error('Detected inlet pores are not coplanar')
    tree = cKDTree(data=inlets)
    Ls = np.unique(np.float64(tree.query(x=outlets)[0]))
    if not np.allclose(Ls, Ls[0]):
        logger.error('A unique value of length could not be found')
    length = Ls[0]
    return length


def reduce_coordination(network, z):
    r"""
    Deletes throats on network to match specified average coordination number

    Parameters
    ----------
    target : Network
        The network whose throats are to be trimmed
    z : scalar
        The desired average coordination number.  It is not possible to specify
        the distribution of the coordination, only the mean value.

    Returns
    -------
    trim : ndarray
        A boolean array with ``True`` values indicating which pores to trim
        (using ``op.topotools.trim``) to obtain the desired average
        coordination number.

    Notes
    -----
    This method first finds the minimum spanning tree of the network using
    random weights on each throat, then assures that these throats are *not*
    deleted, in order to maintain network connectivity.  The list of throats
    to trim is generated randomly from the throats *not* on the spanning tree.

    """
    # Find minimum spanning tree using random weights
    am = network.create_adjacency_matrix(weights=np.random.rand(network.Nt),
                                         triu=False)
    mst = csgraph.minimum_spanning_tree(am, overwrite=True)
    mst = mst.tocoo()

    # Label throats on spanning tree to avoid deleting them
    Ts = network.find_connecting_throat(mst.row, mst.col)
    Ts = np.hstack(Ts)
    network['throat.mst'] = False
    network['throat.mst'][Ts] = True

    # Trim throats not on the spanning tree to acheive desired coordination
    Ts = np.random.permutation(network.throats('mst', mode='nor'))
    del network['throat.mst']
    Ts = Ts[:int(network.Nt - network.Np*(z/2))]
    Ts = network.to_mask(throats=Ts)
    return Ts


def add_reservoir_pore(cls, network, pores, offset=0.1):
    r"""
    Adds a single pore connected to all ``pores`` to act as a reservoir

    This function is mostly needed to make network compatible with the
    Statoil file format, which requires reservoir pores on the inlet and
    outlet faces.

    Parameters
    ----------
    network : Network
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
    # Ps = np.arange(network.Np, network.Np + 1)
    extend(network=network, coords=[new_coord], labels=['reservoir'])
    conns = [[P, network.Np-1] for P in pores]
    # Ts = np.arange(network.Nt, network.Nt + len(conns))
    extend(network=network, conns=conns, labels=['reservoir'])
    # Compute the geometrical properties of the reservoir pore and throats
    # Confirm if network has any geometry props on it
    props = {'throat.length', 'pore.diameter', 'throat.volume'}
    if len(set(network.keys()).intersection(props)) > 0:
        raise Exception('Geometrical properties should be moved to a '
                        + 'geometry object first')
        # or just do this?:  geo = Imported(network=network)
    network.add_model(propname='pore.diameter',
                      model=mods.geometry.pore_size.largest_sphere)
    network.add_model(propname='throat.diameter_temp',
                      model=mods.geometry.throat_size.from_neighbor_pores,
                      mode='min')
    network.add_model(propname='throat.diameter',
                      model=mods.misc.scaled,
                      prop='throat.diameter_temp', factor=0.5)
    network.add_model(propname='throat.volume',
                      model=mods.geometry.throat_volume.cylinder)
    return network.project
