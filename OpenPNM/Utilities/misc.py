import scipy as _sp
import time as _time
from scipy.spatial.distance import cdist as dist

def reflect_pts(coords,nplane):
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

def crop_pts(coords,box):
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
    coords = coords[_sp.any(coords<box[0],axis=1)]
    coords = coords[_sp.any(coords>box[1],axis=1)]
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
    coords = _sp.array(coords,ndmin=1)
    if _sp.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:,0]
    Py = coords[:,1]
    Pz = coords[:,2]

    #Do easy check first, for common coordinate
    if _sp.shape(_sp.unique(Px))[0] == 1:
        return True
    if _sp.shape(_sp.unique(Py))[0] == 1:
        return True
    if _sp.shape(_sp.unique(Pz))[0] == 1:
        return True

    #Perform rigorous check using vector algebra
    n = _sp.array((Px - Px[0],Py - Py[0],Pz - Pz[0])).T
    n0 = _sp.array((Px[-1] - Px[0],Py[-1] - Py[0],Pz[-1] - Pz[0])).T

    n_cross = _sp.cross(n0,n)
    n_dot = _sp.multiply(n0,n_cross)

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
        if quiet == False:
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
            match=False
            for j in output_list:
                if dim == 3:
                    if (i[0]==j[0]) and (i[1]==j[1]) and (i[2]==j[2]):
                        match=True
                elif dim == 2:
                    if (i[0]==j[0]) and (i[1]==j[1]):
                        match=True
                elif dim ==1:
                    if (i[0]==j[0]):
                        match=True
            if match==False:
                output_list.append(i)
    return output_list


class PrintableList(list):
    def __str__(self):
        count = 0
        header = '-'*79
        print('\n')
        print(header)
        self.sort()
        for item in self:
            count = count + 1
            print(count,'\t: ',item)
        return header

class PrintableDict(dict):
    def __str__(self):
        import pprint
        header = '-'*79
        print('\n')
        print(header)
        pprint.pprint(self)
        print(header)
        return ''

class ObjectView(object):
    def __init__(self, d):
        temp = {}
        for item in d:
            if type(d[item][0]) == _sp.bool_:
                key = 'label_'+item.replace('.','_')
            else:
                key = 'prop_'+item.replace('.','_')
            temp[key] =d[item]
        self.__dict__ = temp

def amalgamate_data(objs=[]):
    r"""
    Returns a dictionary containing ALL pore data from all netowrk and/or
    phase objects received as arguments
    """
    if type(objs) is not list:
        objs = list(objs)
    data_amalgamated = {}
    exclusion_list = ['pore.centroid','pore.vertices','throat.centroid','throat.offset_vertices','throat.vertices','throat.normal','throat.perimeter']
    for item in objs:
        mro = [module.__name__ for module in item.__class__.__mro__]
        if 'GenericNetwork' in mro: #if Network object, combine Geometry and Network keys
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
                        dict_name = element+'.'+item.name+'_'+propname
                        if key in ['pore.coords', 'throat.conns', 'pore.all','throat.all']:
                            dict_name = key
                        data_amalgamated.update({dict_name : item[key]})
                except TypeError:
                    #print(key)
                    pass
    return data_amalgamated

def conduit_lengths(network,throats=None,mode='pore'):
    r"""
    Return the respective lengths of the conduit components defined by the throat conns P1 T P2
    mode = 'pore' - uses pore coordinates
    mode = 'centroid' uses pore and throat centroids
    """
    if throats is None:
        throats = network.throats()
    Ps = network['throat.conns']
    pdia = network['pore.diameter']

    if mode ==  'centroid':
        try:
            pcentroids = network['pore.centroid']
            tcentroids = network['throat.centroid']
            plen1 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:,0]]-tcentroids),1))-network['throat.length']/2
            plen2 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:,1]]-tcentroids),1))-network['throat.length']/2
        except KeyError:
            mode = 'pore'
    if mode == 'pore':
        #Find half-lengths of each pore
        pcoords = network['pore.coords']
        #   Find the pore-to-pore distance, minus the throat length
        lengths = _sp.sqrt(_sp.sum(_sp.square(pcoords[Ps[:,0]]-pcoords[Ps[:,1]]),1))-network['throat.length']
        #   Calculate the fraction of that distance from the first pore
        try:
            fractions = pdia[Ps[:,0]]/(pdia[Ps[:,0]]+pdia[Ps[:,1]])
        except:
            fractions = 0.5
        plen1 = lengths*fractions
        plen2 = lengths*(1-fractions)

    return _sp.vstack((plen1,network['throat.length'],plen2)).T[throats]

def clone_object(obj):
    r'''
    Clone an OpenPNM Object

    Parameters
    ----------
    obj : OpenPNM Object
        The object to be cloned can be any OpenPNM Object

    Returns
    -------
    A clone of the specified object is returned, but it retains all its links
    to the objects associated with the original object.  The cloned object is
    not associated with the Network.

    Notes
    -----
    This method is intended to create a disposable object, for instance, to
    receive simulation data without overwriting existing data.

    '''
    cls = obj.__class__.__mro__[0]
    new_obj = cls()
    new_obj.update(obj)
    new_obj._phases = obj._phases
    new_obj._physics = obj._physics
    new_obj._geometries = obj._geometries
    new_obj._net = obj._net
    new_obj._models= obj._models
    return new_obj

def clone_simulation(network,name=None):
    r'''
    Clone an entire Network simulation, including all associated objects.  The
    cloned simulation is numerical identical to, but entirely distinct from the
    original simulation.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network object associated with the simulation to be cloned
    name : string, optional
        A new name can be given to the cloned Network if desired, otherwise it
        will inherit the name of the original Network object.
    '''
    #Clone Network
    cls = network.__class__.__mro__[0]
    new_net = cls(name=name)
    new_net.update(network)
    # Clone associated Geometry
    for item in network._geometries:
        cls = item.__class__.__mro__[0]
        Ps = item.map_pores(pores=item.Ps,target=network)
        Ts = item.map_throats(throats=item.Ts,target=network)
        geom = cls(network=new_net,name=item.name,pores=Ps,throats=Ts)
        geom.update(item)
        geom._models = item._models
    # Clone associated Phases
    for item in network._phases:
        cls = item.__class__.__mro__[0]
        phase = cls(network=new_net,name=item.name)
        phase.update(item)
        phase._models = item._models
    # Repeat Phases to find component phases
    for item in network._phases:
        new_item = new_net.phases(item.name)[0]
        new_item._phases = item._phases
        phase._models = item._models
    # Clone associated Physics
    for item in network._physics:
        cls = item.__class__.__mro__[0]
        phase = item._phases[0]
        Ps = item.map_pores(pores=item.Ps,target=network)
        Ts = item.map_throats(throats=item.Ts,target=network)
        phys = cls(network=new_net,phase=phase,name=item.name,pores=Ps,throats=Ts)
        phys.update(item)
        phys._models = item._models
    return new_net

def _subset(network,pores,name=None):
    r'''
    Create a new sub-network from a given Network, from a list of pores.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network simulation from which a subnet is to be created
    pores : array_like
        A list of pores from which to create the new network
    name : string, optional
        The name to apply to the new network object

    Returns
    -------
    OpenPNM Object
        Returns a new network object

    Notes
    -----
    This is a work in progress

    Examples
    --------
    na
    '''
    import OpenPNM.Utilities.misc as misc
    if name == network.name:
        raise Exception('Subset cannot have same name as parent network')
    # Clone network
    new_net = misc.clone_simulation(network=network,name=name)
    # Add temporary indices to new_net
    new_net['pore.temp_ind'] = network.Ps
    new_net['throat.temp_ind'] = network.Ts
    # Trim cloned network to specific subset
    Ps = ~network.tomask(pores)
    new_net.trim(Ps)
    # Associate cloned network with parent
    new_net._net = network
    # Create labels in parent network
    Ps = new_net['pore.temp_ind']
    network['pore.'+new_net.name] = False
    network['pore.'+new_net.name][Ps] = True
    Ts = new_net['throat.temp_ind']
    network['throat.'+new_net.name] = False
    network['throat.'+new_net.name][Ts] = True

    return new_net
