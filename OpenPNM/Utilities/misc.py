import scipy as _sp

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