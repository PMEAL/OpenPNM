"""
module __Template__: Generate cubic networks from domain templates
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM
import matplotlib as plt
import scipy as sp
import scipy.ndimage as spim
import numpy as np
from OpenPNM.Network.__Template__ import Template

class Sphere(Template):
    r"""
    This class contains the methods to create a cubic network with an spherical 
    domain of specified radius.
    
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
    """

    def __init__(self, **kwargs):
        super(Sphere,self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__,": ","Execute constructor")
    
    def generate(self, **params):
        r'''
        radius: int
            The physical radius of the netowrk.  The number of pores determined
            from this and the lattice_spacing parameter.
        lattice_spacing : float
            The lattice constant for the network, used to scale distance between pores.
        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        #create spherical template, then call normal template procedure
        self._radius = params['radius']
        Nr = sp.around(self._radius/params['lattice_spacing'])
        temp = sp.ones((2*Nr, 2*Nr, 2*Nr))
        temp[Nr, Nr, Nr] = 0
        temp = spim.distance_transform_edt(temp)
        temp = temp < Nr
        params['template'] = temp
        #Call standard generation protocol
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        self._add_boundaries()
        self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")        
        return self
        

if __name__ == '__main__':
    pn = OpenPNM.Network.Sphere(name='sphere_1').generate(radius=5, lattice_spacing=1)
    print(pn.name)







