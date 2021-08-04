import numpy as np
import scipy as sp
import scipy.ndimage as spim
from scipy.sparse import csgraph
from scipy.spatial import ConvexHull
from openpnm.utils import logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


def isoutside(coords, shape):
    r"""
    Identifies points that lie outside the specified shape

    Parameters
    ----------
    coords : array_like
        The coordinates which are to be checked
    shape : array_like
        The shape of the domain beyond which points should be trimmed.
        The argument is treated as follows:

        **sphere** : If a scalar or single element list is received, it's
        treated as the radius [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as
        the radius and height of a cylinder [r, z] whose central axis
        starts at [0, 0, 0] and extends in the positive z-direction.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner
        lies at [0, 0, 0].

    Returns
    -------
    An Np-long mask of ``True`` values indicating pores that lie outside the
    domain.

    """
    # Label external pores for trimming below
    if len(shape) == 1:  # Spherical
        # Find external points
        r = np.sqrt(np.sum(coords**2, axis=1))
        Ps = r > shape[0]
    elif len(shape) == 2:  # Cylindrical
        # Find external pores outside radius
        r = np.sqrt(np.sum(coords[:, [0, 1]]**2, axis=1))
        Ps = r > shape[0]
        # Find external pores above and below cylinder
        if shape[1] > 0:
            Ps = Ps + (coords[:, 2] > shape[1])
            Ps = Ps + (coords[:, 2] < 0)
        else:
            pass
    elif len(shape) == 3:  # Rectilinear
        shape = np.array(shape, dtype=float)
        try:
            lo_lim = shape[:, 0]
            hi_lim = shape[:, 1]
        except IndexError:
            lo_lim = np.array([0, 0, 0])
            hi_lim = shape
        Ps1 = np.any(coords > hi_lim, axis=1)
        Ps2 = np.any(coords < lo_lim, axis=1)
        Ps = Ps1 + Ps2
    return Ps


def rotate_coords(network, a=0, b=0, c=0, R=None):
    r"""
    Rotates coordinates a given amount about each axis

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose pore coordinates should be transformed
    a : scalar
        The amount in degrees to rotate about the x-axis
    b : scalar
        The amount in degrees to rotate about the y-axis
    c : scalar
        The amount in degrees to rotate about the z-axis
    R : array_like
        Rotation matrix.  Must be a 3-by-3 matrix since pore coordinates are
        always in 3D.  If this is given then the other individual arguments
        are ignored.

    See Also
    --------
    rotate_coords

    Notes
    -----
    It is possible to rotate about any of the three axes by specifying ``a``,
    ``b``, and/or ``c``.  In this case each rotation is applied in sequence.

    """
    if R is None:
        if a:
            R = np.array([[1, 0, 0],
                          [0, np.cos(np.deg2rad(a)), -np.sin(np.deg2rad(a))],
                          [0, np.sin(np.deg2rad(a)), np.cos(np.deg2rad(a))]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
        if b:
            R = np.array([[np.cos(np.deg2rad(b)), 0, -np.sin(np.deg2rad(b))],
                          [0, 1, 0],
                          [np.sin(np.deg2rad(b)), 0, np.cos(np.deg2rad(b))]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
        if c:
            R = np.array([[np.cos(np.deg2rad(c)), -np.sin(np.deg2rad(c)), 0],
                          [np.sin(np.deg2rad(c)), np.cos(np.deg2rad(c)), 0],
                          [0, 0, 1]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
    else:
        network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                              axes=(1, 1))


def shear_coords(network, ay=0, az=0, bx=0, bz=0, cx=0, cy=0, S=None):
    r"""
    Shears the coordinates a given amount about along axis

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose pore coordinates should be transformed
    ay : scalar
        The factor by which to shear along the x-axis as a function of y
    az : scalar
        The factor by which to shear along the x-axis as a function of z
    bx : scalar
        The factor by which to shear along the y-axis as a function of x
    bz : scalar
        The factor by which to shear along the y-axis as a function of z
    cx : scalar
        The factor by which to shear along the z-axis  as a function of x
    cy : scalar
        The factor by which to shear along the z-axis as a function of y
    S : array_like
        The shear matrix.  Must be a 3-by-3 matrix since pore coordinates are
        always in 3D.  If this is given then the other individual arguments
        are ignored.

    See Also
    --------
    rotate_coords

    Notes
    -----
    The shear along the i *th* -axis is given as i\* = i + aj.  This means
    the new i coordinate is the old one plus some linear factor *a* in the
    j *th* direction.

    The values of ``a``, ``b``, and ``c`` are essentially the inverse of the
    slope to be formed by the neighboring layers of sheared pores.  A value of
    0 means no shear, and neighboring points are stacked directly on top of
    each other; a value of 1 means they form a 45 degree diagonal, and so on.

    If ``S`` is given, then is should be of the form:

    ::

        S = [[1 , ay, az],
             [bx, 1 , bz],
             [cx, cy, 1 ]]

        where any of the off-diagonal components can be 0

    """
    if S is None:
        S = np.array([[1, ay, az],
                      [bx, 1, bz],
                      [cx, cy, 1]])
    network['pore.coords'] = (S@network['pore.coords'].T).T


def trim(network, pores=[], throats=[]):
    '''
    Remove pores or throats from the network

    Parameters
    ----------
    network : OpenPNM Network Object
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

    '''
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
            logger.info('Removing ALL throats from network')
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
        else:
            Ps = obj.map_pores(pores=Pkeep, origin=network)
            Ts = obj.map_throats(throats=Tkeep, origin=network)
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
    r'''
    Add pores or throats to the network from a list of coords or conns.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network to which pores or throats should be added
    coords : array_like
        The coordinates of the pores to add.  These will be appended to the
        'pore.coords' array so should be of shape N-by-3, where N is the
        number of pores in the list.
    conns : array_like
        The throat connections to add.  These will be appended to the
        'throat.conns' array so should be of shape N-by-2.  Note that the
        numbering must point to existing pores.
    labels : string, or list of strings, optional
        A list of labels to apply to the new pores and throats

    '''
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


def reduce_coordination(network, z):
    r"""
    Deletes throats on network to match specified average coordination number

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose throats are to be trimmed
    z : scalar
        The desire average coordination number.  It is not possible to specify
        the distribution of the coordination, only the mean value.

    Notes
    -----
    This method first finds the minimum spanning tree of the network using
    random weights on each throat, then assures that these throats are *not*
    deleted, in order to maintain network connectivity.

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
    Ts = Ts[:int(network.Nt - network.Np*(z/2))]
    trim(network=network, throats=Ts)


def label_faces(network, tol=0.0, label='surface'):
    r"""
    Finds pores on the surface of the network and labels them according to
    whether they are on the *top*, *bottom*, etc.  This function assumes the
    network is cubic in shape (i.e. with six flat sides)

    Parameters
    ----------
    network : OpenPNM Network object
        The network to apply the labels

    tol : scalar
        The tolerance for defining what counts as a surface pore, which is
        specifically meant for random networks.  All pores with ``tol`` of
        the maximum or minimum along each axis are counts as pores.  The
        default is 0.

    label : string
        An identifying label to isolate the pores on the faces of the network.
        The default is 'surface'.  Surface pores can be found using
        ``find_surface_pores``.

    """
    label = label.split('.', 1)[-1]
    if 'pore.'+label not in network.labels():
        find_surface_pores(network, label=label)
    Psurf = network['pore.'+label]
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
        network['pore.right'] = (crds[:, 0] >= (xmax - tol*xspan)) * Psurf
    if dims[1]:
        network['pore.front'] = (crds[:, 1] <= (ymin + tol*yspan)) * Psurf
        network['pore.back'] = (crds[:, 1] >= (ymax - tol*yspan)) * Psurf
    if dims[2]:
        network['pore.top'] = (crds[:, 2] >= (zmax - tol*zspan)) * Psurf
        network['pore.bottom'] = (crds[:, 2] <= (zmin + tol*zspan)) * Psurf


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
        3 x N array of the marker coordinates to use in the triangulation.  The
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


def dimensionality(network):
    r"""
    Checks the dimensionality of the network

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose dimensionality is to be checked

    Returns
    -------
    dims : list
        A  3-by-1 array containing ``True`` for each axis that contains
        multiple values, indicating that the pores are spatially distributed
        in that direction.

    """
    xyz = network["pore.coords"]
    eps = np.finfo(float).resolution
    dims_unique = [not np.allclose(xk, xk.mean(), atol=0, rtol=eps) for xk in xyz.T]
    return np.array(dims_unique)


def clone_pores(network, pores, labels=['clone'], mode='parents'):
    r"""
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
    network : OpenPNM Network Object
        The network to which all the other networks should be added.

    donor : OpenPNM Network Object or list of Objects
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
                        logger.debug('Adding ' + key + ' to donor')
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
    r'''
    Stitches a second a network to the current network.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network to which to donor Network will be attached

    donor : OpenPNM Network Object
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

    label_stitches : str or list of strings
        The label to apply to the newly created 'stitch' throats.  The
        defaul is 'stitched'.  If performing multiple stitches in a row it
        might be helpful to the throats created during each step uniquely
        for later identification.

    len_max : float
        Set a length limit on length of new throats

    method : string (default = 'nearest')
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

    '''
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

    if len(network.project.geometries()) > 0:
        logger.warning(str(conns.shape[0]) + ' newly created throats are not '
                       + 'assigned to a geometry')

    # Remove donor from Workspace, if present
    # This check allows for the reuse of a donor Network multiple times
    for sim in list(ws.values()):
        if donor in sim:
            del ws[sim.name]


def stitch_pores(network, pores1, pores2, mode='gabriel'):
    r"""
    Stitches together pores in a network with disconnected clusters

    Parameter
    ---------
    network : OpenPNM Network
        The network to operate upon
    pores1 and pores2: array_like
        The pore indices of the disconnected clusters to be joined
    mode : str
        Dictates which tesselation method is used to identify which pores to
        stitch together.  Options are 'gabriel' (default) or 'delaunay'.

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
    r'''
    Returns the possible connections between two groups of pores, and optionally
    makes the connections.

    See ``Notes`` for advanced usage.

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

    '''
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


def find_pore_to_pore_distance(network, pores1=None, pores2=None):
    r'''
    Find the distance between all pores on set 1 to each pore in set 2

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
    dist : array_like
        A distance matrix with ``len(pores1)`` rows and ``len(pores2)`` columns.
        The distance between pore *i* in ``pores1`` and *j* in ``pores2`` is
        located at *(i, j)* and *(j, i)* in the distance matrix.

    Notes
    -----
    This function computes and returns a distance matrix, which is a dense
    matrix of size Np_1 by Np_2, so can get large.  For distances between
    larger sets a KD-tree approach would be better, which is available in
    ``scipy.spatial``.

    '''
    from scipy.spatial.distance import cdist
    p1 = np.array(pores1, ndmin=1)
    p2 = np.array(pores2, ndmin=1)
    coords = network['pore.coords']
    return cdist(coords[p1], coords[p2])


def filter_pores_by_z(network, pores, z=1):
    r"""
    Find pores with a given number of neighbors

    Parameters
    ----------
    network : OpenPNM Network object
        The network on which the query is to be performed
    pores : array_like
        The pores to be filtered
    z : int
        The coordination number to filter by

    Returns
    -------
    pores : array_like
        The pores which have the specified coordination number

    """
    pores = network._parse_indices(pores)
    Nz = network.num_neighbors(pores=pores)
    orphans = np.where(Nz == z)[0]
    hits = pores[orphans]
    return hits


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

    '''
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
        raise Exception('The network has subdivided pores, so the method \
                         does not support another subdivision')
    # Assigning right shape and division
    if np.size(shape) != 2 and np.size(shape) != 3:
        raise Exception('Subdivide not implemented for Networks other than 2D and 3D')
    if np.size(shape) == 3 and 1 not in shape:
        div = np.array(shape, ndmin=1)
        single_dim = None
    else:
        single_dim = np.where(np.array(network.shape) == 1)[0]
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
    networkspacing = network.spacing
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
    if np.sum(occluded_ts) > 0:
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
    im : array_like
        A Numpy array containing 1's to demarcate the desired shape, and 0's
        elsewhere.

    """
    rmax = np.array(outer_radius, ndmin=1)
    rmin = np.array(inner_radius, ndmin=1)
    ind = 2 * rmax - 1
    coord = np.indices((ind * np.ones(dim, dtype=int)))
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


def template_sphere_shell(outer_radius, inner_radius=0, dim=3):
    r"""
    This method generates an image array of a sphere-shell.

    It is useful for passing to Cubic networks as a ``template`` to make
    spherical shaped networks.

    Parameters
    ----------
    outer_radius : int
        Number of nodes in the outer radius of the sphere.

    inner_radius : int
        Number of nodes in the inner radius of the shell.  a value of 0 will
        result in a solid sphere.

    dim : scalar
        Controls the number of dimensions of the result.  3 returns a sphere,
        while 2 returns a disk.

    Returns
    -------
    im : array_like
        A Numpy array containing 1's to demarcate the sphere-shell, and 0's
        elsewhere.

    """
    img = _template_sphere_disc(dim=dim, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return img


def template_cylinder_annulus(height, outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a disc-ring.

    It is useful for passing to Cubic networks as a ``template`` to make
    circular-shaped 2D networks.

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
    im : array_like
        A Numpy array containing 1's to demarcate the disc-ring, and 0's
        elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    img = np.tile(np.atleast_3d(img), reps=height)
    return img


def generate_base_points(num_points, domain_size, density_map=None,
                         reflect=True):
    r"""
    Generates a set of base points for passing into the Tessellation-based
    Network classes.  The points can be distributed in spherical, cylindrical,
    or rectilinear patterns, as well as 2D and 3D (disks and squares).

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
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.

    density_map : array, optional
        A an array that contains fractional values (0 < i < 1) indicating the
        liklihood that a point in that region should be kept.  The size of this
        array can be anything, but the shape must match the ``domain_size``;
        that is for a 3D network the shape of the ``density_map`` can be
        [10, 10, 10] or [50, 50, 50], depending on how important the resolution
        of the density distribution is.  For a 2D network the ``density_map``
        should be [10, 10].

        When specifying a custom probabiliy map is it recommended to also set
        values outside the given domain to zero.  If not, then the correct
        shape will still be returned, but with too few points in it.

    reflect : boolean
        If True, the the base points are generated as specified, the reflected
        about each face of the domain.  This essentially tricks the
        tessellation functions into creating smooth flat faces at the
        boundaries once these excess pores are trimmed.

    Notes
    -----
    The reflection approach tends to create larger pores near the surfaces, so
    it might be necessary to use the ``density_map`` argument to specify a
    slightly higher density of points near the surfaces.

    The ``Voronoi``, ``Delaunay``, ``Gabriel``, and ``DelunayVoronoiDual``
    classes can *techncially* handle base points with spherical or cylindrical
    domains, but the reflection across round surfaces does not create perfect
    Voronoi cells so the surfaces will not be smooth.


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
    >>> im = np.ones([21, 21, 21], dtype=int)
    >>> im[10, 10, 10] = 0
    >>> im = spim.distance_transform_edt(im) <= 20  # Create sphere of 1's
    >>> prob = spim.distance_transform_edt(im)
    >>> prob = prob / np.amax(prob)  # Normalize between 0 and 1
    >>> pts = op.topotools.generate_base_points(num_points=50,
    ...                                         domain_size=[1, 1, 1],
    ...                                         density_map=prob)
    >>> net = op.network.DelaunayVoronoiDual(points=pts, shape=[1, 1, 1])

    """
    def _try_points(num_points, prob):
        prob = np.atleast_3d(prob)
        prob = np.array(prob)/np.amax(prob)  # Ensure prob is normalized
        base_pts = []
        N = 0
        while N < num_points:
            pt = np.random.rand(3)  # Generate a point
            # Test whether to keep it or not
            [indx, indy, indz] = np.floor(pt*np.shape(prob)).astype(int)
            if np.random.rand(1) <= prob[indx][indy][indz]:
                base_pts.append(pt)
                N += 1
        base_pts = np.array(base_pts)
        return base_pts

    if len(domain_size) == 1:  # Spherical
        domain_size = np.array(domain_size)
        r = domain_size[0]
        if density_map is None:
            # Make an image of a sphere filled with ones and use _try_points
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, 20] = 0
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to spherical coordinates
        X, Y, Z = np.array(base_pts - [0.5, 0.5, 0.5]).T
        r = 2*np.sqrt(X**2 + Y**2 + Z**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        phi = 2*np.arctan(np.sqrt(X**2 + Y**2)/Z)
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, phi] = [r[inds], theta[inds], phi[inds]]
        # Reflect base points across perimeter
        if reflect:
            r, theta, phi = reflect_base_points(np.vstack((r, theta, phi)),
                                                domain_size)
        # Convert to Cartesean coordinates
        X = r*np.cos(theta)*np.sin(phi)
        Y = r*np.sin(theta)*np.sin(phi)
        Z = r*np.cos(phi)
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 2:  # Cylindrical or Disk
        domain_size = np.array(domain_size)
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, :] = 0
            if domain_size[1] == 0:  # Disk
                density_map = density_map[:, :, 0]
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to cylindrical coordinates
        X, Y, Z = np.array(base_pts - [0.5, 0.5, 0]).T  # Center on z-axis
        r = 2*np.sqrt(X**2 + Y**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        z = Z*domain_size[1]
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        inds = ~((z > domain_size[1]) + (z < 0))
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        if reflect:
            r, theta, z = reflect_base_points(np.vstack([r, theta, z]),
                                              domain_size)
        # Convert to Cartesean coordinates
        X = r*np.cos(theta)
        Y = r*np.sin(theta)
        Z = z
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 3:  # Cube or square
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            if domain_size[2] == 0:
                density_map = density_map[:, :, 0]
        base_pts = _try_points(num_points, density_map)
        base_pts = base_pts*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts


def reflect_base_points(base_pts, domain_size):
    r'''
    Helper function for relecting a set of points about the faces of a
    given domain.

    Parameters
    ----------
    base_pts : array_like
        The coordinates of the base_pts to be reflected in the coordinate
        system corresponding to the the domain as follows:

        **spherical** : [r, theta, phi]
        **cylindrical** or **circular** : [r, theta, z]
        **rectangular** or **square** : [x, y, z]

    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.

    '''
    domain_size = np.array(domain_size)
    if len(domain_size) == 1:
        r, theta, phi = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        phi = np.hstack([phi, phi])
        base_pts = np.vstack((r, theta, phi))
    if len(domain_size) == 2:
        r, theta, z = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        z = np.hstack([z, z])
        if domain_size[1] != 0:  # If not a disk
            r = np.hstack([r, r, r])
            theta = np.hstack([theta, theta, theta])
            z = np.hstack([z, -z, 2*domain_size[1]-z])
        base_pts = np.vstack((r, theta, z))
    elif len(domain_size) == 3:
        Nx, Ny, Nz = domain_size
        # Reflect base points about all 6 faces
        orig_pts = base_pts
        base_pts = np.vstack((base_pts,
                              [-1, 1, 1] * orig_pts + [2.0 * Nx, 0, 0]))
        base_pts = np.vstack((base_pts, [-1, 1, 1] * orig_pts))
        base_pts = np.vstack((base_pts,
                              [1, -1, 1] * orig_pts + [0, 2.0 * Ny, 0]))
        base_pts = np.vstack((base_pts, [1, -1, 1] * orig_pts))
        if domain_size[2] != 0:
            base_pts = np.vstack((base_pts,
                                  [1, 1, -1] * orig_pts + [0, 0, 2.0 * Nz]))
            base_pts = np.vstack((base_pts, [1, 1, -1] * orig_pts))
    return base_pts


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
    Ps = np.array(pores, ndmin=1)
    if Ps.dtype is bool:
        Ps = network.toindices(Ps)
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
    results : bool
        A boolean value of whether given points are coplanar (``True``) or
        not (``False``)
    '''
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
    network : GenericNetwork
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
