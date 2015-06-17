import scipy as _sp
import time as _time
from scipy.spatial.distance import cdist as dist


def reflect_pts(coords, nplane):
    r'''
    Reflect points across the given plane

    Parameters
    ----------
    coords : array_like
        An Np x ndims array off [x,y,z] coordinates

    nplane : array_like
        A vector of length ndims, specifying the normal to the plane.  The tail
        of the vector is assume to lie on the plane, and the reflection will
        be applied in the direction of the vector.

    Returns
    -------
    coords : array_like
        An Np x ndmins vector of reflected point, not including the input points

    '''
    pass


def crop_pts(coords, box):
    r'''
    Drop all points lying outside the box

    Parameters
    ----------
    coords : array_like
        An Np x ndims array off [x,y,z] coordinates

    box : array_like
        A 2 x ndims array of diametrically opposed corner coordintes

    Returns
    -------
    coords : array_like
        Inputs coordinates with outliers removed

    Notes
    -----
    This needs to be made more general so that an arbitray cuboid with any
    orientation can be supplied, using Np x 8 points
    '''
    coords = coords[_sp.any(coords < box[0], axis=1)]
    coords = coords[_sp.any(coords > box[1], axis=1)]
    return coords


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
    A boolean value of whether given points are colplanar (True) or not (False)
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
    n = _sp.array((Px - Px[0], Py - Py[0], Pz - Pz[0])).T
    n0 = _sp.array((Px[-1] - Px[0], Py[-1] - Py[0], Pz[-1] - Pz[0])).T

    n_cross = _sp.cross(n0, n)
    n_dot = _sp.multiply(n0, n_cross)

    if _sp.sum(_sp.absolute(n_dot)) == 0:
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


def amalgamate_data(objs=[]):
    r"""
    Returns a dictionary containing ALL pore data from all netowrk and/or
    phase objects received as arguments
    """
    if type(objs) is not list:
        objs = list(objs)
    data_amalgamated = {}
    exclusion_list = ['pore.centroid', 'pore.vertices', 'throat.centroid',
                      'throat.offset_vertices', 'throat.vertices', 'throat.normal',
                      'throat.perimeter', 'pore.vert_index', 'throat.vert_index']
    for item in objs:
        mro = [module.__name__ for module in item.__class__.__mro__]
        # If Network object, combine Geometry and Network keys
        if 'GenericNetwork' in mro:
            keys = []
            for key in item.keys():
                keys.append(key)
            for geom in item._geometries:
                for key in geom.keys():
                    if key not in keys:
                        keys.append(key)
        else:
            if 'GenericPhase' in mro:
                keys = []
                for key in item.keys():
                    keys.append(key)
                for physics in item._physics:
                    for key in physics.keys():
                        if key not in keys:
                            keys.append(key)
        keys.sort()
        for key in keys:
            if key not in exclusion_list:
                try:
                    if _sp.amax(item[key]) < _sp.inf:
                        element = key.split('.')[0]
                        propname = key.split('.')[1]
                        dict_name = element + '.' + item.name + '_' + propname
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
        #   Calculate the fraction of that distance from the first pore
        try:
            fractions = pdia[Ps[:, 0]]/(pdia[Ps[:, 0]] + pdia[Ps[:, 1]])
        except:
            fractions = 0.5
        plen1 = lengths*fractions
        plen2 = lengths*(1-fractions)

    return _sp.vstack((plen1, network['throat.length'], plen2)).T[throats]
