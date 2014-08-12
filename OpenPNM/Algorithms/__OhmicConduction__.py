"""

module __OhmicConduction__:
===============================================================================

"""

import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class OhmicConduction(GenericLinearTransport):
    r'''

    '''

    def __init__(self,**kwargs):
        r'''
        '''
        super(OhmicConduction,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,fluid,conductance='electronic_conductance',quantity='voltage',**params):
        r'''
        '''  
        self._logger.info("Setup "+self.__class__.__name__)        
        super(OhmicConduction,self).setup(fluid=fluid,conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()