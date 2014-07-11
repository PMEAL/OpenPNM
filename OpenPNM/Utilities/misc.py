import scipy as _sp
import time as _time
from scipy.spatial.distance import cdist as dist


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
        header = '-'*50
        print(header)
        for item in self:
            count = count + 1
            print(count,'\t: ',item)
        return header























        
