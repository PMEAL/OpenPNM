"""
===============================================================================
DelaunayCubic: Generate random networks based on Delaunay Tessellations
===============================================================================

"""
import OpenPNM
import scipy as sp
import sys
import numpy as np
from OpenPNM.Network.__Delaunay__ import Delaunay


class DelaunayCubic(Delaunay):
    r"""
    This class contains the methods for creating a *Delaunay* network topology
    based connecting pores with a Delaunay tessellation.
    
    The Subclass Generates Points on a cubic lattice and then perturbs them to prevent degeneracy
    
    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network
        
    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
        
    loggername : string
        Overwrite the name of the logger, which defaults to the class name

    Examples
    --------
    >>> pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[0.0001,0.0001,0.0001],name='net')
    >>> pn.num_pores()
    100

    """

    def __init__(self,shape=None, template=None, spacing=[1,1,1],jiggle_factor=0.1,**kwargs):
        '''
        Create Delauny network object
        '''
        if shape != None:
            self._arr = np.atleast_3d(np.empty(shape))
        elif template != None:
            self._arr = sp.array(template,ndmin=3,dtype=bool)
        
        self._shape = sp.shape(self._arr)  # Store original network shape
        self._spacing = sp.asarray(spacing)  # Store network spacing instead of calculating it
        self._num_pores = np.prod(np.asarray(self._shape))
        self._domain_size = np.asarray(self._shape)*self._spacing       
        self._jiggle_factor = jiggle_factor
        
        super(DelaunayCubic,self).__init__(num_pores=self._num_pores,domain_size=self._domain_size,**kwargs)
        
    
    def _generate_pores(self):
        r"""
        Generate the pores with numbering scheme.
        """
        
        points = np.array([i for i,v in np.ndenumerate(self._arr)], dtype=float)
        points += 0.5
        
        "-----------------------------------------------------------------------------------"        
        "2D Orthorhombic adjustment - shift even rows back a bit and odd rows forward a bit"
        "   0   0   0   "
        " 0   0   0   0 "
        "   0   0   0   "
        #shift_y=np.array([0,0.25,0])
        #shift_x=np.array([0.25,0,0])
        #points[(points[:,0] % 2 == 0)] -= shift_y
        #points[(points[:,2] % 2 != 0)] -= shift_x
        #points[(points[:,0] % 2 != 0)] += shift_y
        #points[(points[:,2] % 2 == 0)] += shift_x
        "-----------------------------------------------------------------------------------" 
        body_points=[]
        for i in range(1,self._shape[0]):
            for j in range(1,self._shape[1]):
                for k in range(1,self._shape[2]):
                    body_points.append([i,j,k])
        body_points = np.asarray(body_points)
        points = np.concatenate((points,body_points))
        
        jiggle = (np.random.rand(len(points),3)-0.5)*self._jiggle_factor
        points += jiggle
        points *= self._spacing
        
        self['pore.coords']  = points
        self._logger.debug(sys._getframe().f_code.co_name+": End of method") 
           


if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
