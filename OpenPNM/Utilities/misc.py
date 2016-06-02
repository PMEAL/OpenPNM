import scipy as _sp
import time as _time
import scipy.sparse as _sprs
import OpenPNM as _op
from scipy.spatial.distance import cdist as dist


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
    >>> import OpenPNM
    >>> import OpenPNM.Utilities.misc as misc
    >>> pn = OpenPNM.Network.Cubic(shape=[3, 3, 3])
    >>> a = misc.find_path(network=pn, pore_pairs=[[0, 4], [0, 10]])
    >>> a['pores']
    [array([0, 1, 4]), array([ 0,  1, 10])]
    >>> a['throats']
    [array([ 0, 19]), array([ 0, 37])]
    """
    Ps = _sp.array(pore_pairs, ndmin=2)
    if weights is None:
        weights = _sp.ones_like(network.Ts)
    graph = network.create_adjacency_matrix(data=weights,
                                            sprsfmt='csr',
                                            dropzeros=False)
    paths = _sprs.csgraph.dijkstra(csgraph=graph,
                                   indices=Ps[:, 0],
                                   return_predecessors=True)[1]
    pores = []
    throats = []
    for row in range(0, _sp.shape(Ps)[0]):
        j = Ps[row][1]
        ans = []
        while paths[row][j] > -9999:
            ans.append(j)
            j = paths[row][j]
        ans.append(Ps[row][0])
        ans.reverse()
        pores.append(_sp.array(ans))
        throats.append(network.find_neighbor_throats(pores=ans,
                                                     mode='intersection'))
    pdict = _op.Base.Tools.PrintableDict
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
    coords = _sp.array(coords, ndmin=1)
    if _sp.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:, 0]
    Py = coords[:, 1]
    Pz = coords[:, 2]

    # Do easy check first, for common coordinate
    if _sp.shape(_sp.unique(Px))[0] == 1:
        return True
    if _sp.shape(_sp.unique(Py))[0] == 1:
        return True
    if _sp.shape(_sp.unique(Pz))[0] == 1:
        return True

    # Perform rigorous check using vector algebra
    n1 = _sp.array((Px[1] - Px[0], Py[1] - Py[0], Pz[1] - Pz[0])).T
    n2 = _sp.array((Px[2] - Px[1], Py[2] - Py[1], Pz[2] - Pz[1])).T
    n = _sp.cross(n1, n2)
    r = _sp.array((Px[1:-1] - Px[0], Py[1:-1] - Py[0], Pz[1:-1] - Pz[0]))

    n_dot = _sp.dot(n, r)

    if _sp.sum(n_dot) == 0:
        return True
    else:
        return False


def tic():
    r'''
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.
    '''
    global _startTime_for_tictoc
    _startTime_for_tictoc = _time.time()


def toc(quiet=False):
    r'''
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    Parameters
    ----------
    quiet : Boolean
        If False (default) then a message is output to the console.  If True
        the message is not displayed and the elapsed time is returned.
    '''
    if '_startTime_for_tictoc' in globals():
        t = _time.time() - _startTime_for_tictoc
        if quiet is False:
            print('Elapsed time in seconds: ', t)
        else:
            return t
    else:
        print("Toc: start time not set")


def unique_list(input_list):
    r"""
    For a given list (of points) remove any duplicates
    """
    output_list = []
    if len(input_list) > 0:
        dim = _sp.shape(input_list)[1]
        for i in input_list:
            match = False
            for j in output_list:
                if dim == 3:
                    if i[0] == j[0] and i[1] == j[1] and i[2] == j[2]:
                        match = True
                elif dim == 2:
                    if i[0] == j[0] and i[1] == j[1]:
                        match = True
                elif dim == 1:
                    if i[0] == j[0]:
                        match = True
            if match is False:
                output_list.append(i)
    return output_list


def amalgamate_data(objs=[], delimiter='_'):
    r"""
    Returns a dictionary containing ALL pore data from all netowrk and/or
    phase objects received as arguments

    Parameters
    ----------
    obj : list of OpenPNM objects
        The network and Phase objects whose data should be amalgamated into a
        single dict

    delimiter : string
        The delimiter to place between the prop name and the object name.  For
        instance \'pore.air_molar_density\' or \'pore.air|molar_density'\.  The
        use of underscores can be problematic for reloading the data since they
        are also used in multiple word properties.  The default is '_' for
        backwards compatibility, but the '|' option is preferred.

    Returns
    -------
    A standard Python dict containing all the data from the supplied OpenPNM
    objects
    """
    if type(objs) is not list:
        objs = list(objs)
    data_amalgamated = {}
    dlim = delimiter
    exclusion_list = ['pore.centroid', 'pore.vertices', 'throat.centroid',
                      'throat.offset_vertices', 'throat.vertices', 'throat.normal',
                      'throat.perimeter', 'pore.vert_index', 'throat.vert_index']
    for item in objs:
        mro = [module.__name__ for module in item.__class__.__mro__]
        # If Network object, combine Geometry and Network keys
        if 'GenericNetwork' in mro:
            keys = []
            for key in list(item.keys()):
                keys.append(key)
            for geom in item._geometries:
                for key in list(geom.keys()):
                    if key not in keys:
                        keys.append(key)
        else:
            if 'GenericPhase' in mro:
                keys = []
                for key in list(item.keys()):
                    keys.append(key)
                for physics in item._physics:
                    for key in list(physics.keys()):
                        if key not in keys:
                            keys.append(key)
        keys.sort()
        for key in keys:
            if key not in exclusion_list:
                try:
                    if _sp.amax(item[key]) < _sp.inf:
                        element = key.split('.')[0]
                        propname = key.split('.')[1]
                        dict_name = element + '.' + item.name + dlim + propname
                        if key in ['pore.coords', 'throat.conns',
                                   'pore.all', 'throat.all']:
                            dict_name = key
                        data_amalgamated.update({dict_name: item[key]})
                except TypeError:
                    pass
    return data_amalgamated


def conduit_lengths(network, throats=None, mode='pore'):
    r"""
    Return the respective lengths of the conduit components defined by the throat
    conns P1 T P2
    mode = 'pore' - uses pore coordinates
    mode = 'centroid' uses pore and throat centroids
    """
    if throats is None:
        throats = network.throats()
    Ps = network['throat.conns']
    pdia = network['pore.diameter']

    if mode == 'centroid':
        try:
            pcentroids = network['pore.centroid']
            tcentroids = network['throat.centroid']
            if _sp.sum(_sp.isnan(pcentroids)) + _sp.sum(_sp.isnan(tcentroids)) > 0:
                mode = 'pore'
            else:
                plen1 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:, 0]] -
                                         tcentroids), 1))-network['throat.length']/2
                plen2 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:, 1]] -
                                         tcentroids), 1))-network['throat.length']/2
        except KeyError:
            mode = 'pore'
    if mode == 'pore':
        # Find half-lengths of each pore
        pcoords = network['pore.coords']
        # Find the pore-to-pore distance, minus the throat length
        lengths = _sp.sqrt(_sp.sum(_sp.square(pcoords[Ps[:, 0]] -
                                   pcoords[Ps[:, 1]]), 1)) - network['throat.length']
        lengths[lengths < 0.0] = 2e-9
        # Calculate the fraction of that distance from the first pore
        try:
            fractions = pdia[Ps[:, 0]]/(pdia[Ps[:, 0]] + pdia[Ps[:, 1]])
            # Don't allow zero lengths
#            fractions[fractions == 0.0] = 0.5
#            fractions[fractions == 1.0] = 0.5
        except:
            fractions = 0.5
        plen1 = lengths*fractions
        plen2 = lengths*(1-fractions)

    return _sp.vstack((plen1, network['throat.length'], plen2)).T[throats]
