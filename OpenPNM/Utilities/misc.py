import scipy as _sp
import time as _time
from scipy.spatial.distance import cdist as dist
import OpenPNM.Utilities.transformations as tr

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
    #angles=[]
    #for vec in n_cross[2:len(n_cross)]:
    #    angles.append(180*(tr.angle_between_vectors(n_cross[1],vec,directed=False)/_sp.pi))
    #angles=_sp.asarray(angles)
    #if angles.mean() < 20:
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
    if type(objs) != list:
        objs = list(objs)
    data_amalgamated = {}
    exclusion_list = ['pore.centroid','pore.vertices','throat.centroid','throat.offset_verts','throat.verts','throat.normals','throat.perimeter']
    for item in objs:
        if item.__module__.split('.')[1] == 'Network': #if network object, combine geometry and network keys
            keys = []
            for key in item.keys():
                keys.append(key)
            for geom in item._geometries:
                for key in geom.keys():
                    if key not in keys:
                        keys.append(key)
        else:
            if item.__module__.split('.')[1] == 'Phases':
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
                if _sp.amax(item[key]) < _sp.inf:
                    dict_name = item.name+'.'+key
                    data_amalgamated.update({dict_name : item[key]})
    return data_amalgamated




















        
