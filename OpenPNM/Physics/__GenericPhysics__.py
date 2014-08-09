"""
module Physics
===============================================================================

"""
import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import scipy as sp

class GenericPhysics(OpenPNM.Core):
    r"""
    Generic class to generate Physics objects  

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Physics should be attached
        
    fluid : OpenPNM Fluid object 
        The Fluid object to which this Physics applies
    
    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.
    
    name : str, optional
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.  
    
    """

    def __init__(self,network,fluid,pores=[],throats=[],name=None,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        #Append objects self for internal access
        self._net = network
        self.name = name
        
        #Append self to other objects
        network._physics.append(self)
        fluid._physics.append(self)
        self._fluids.append(fluid)
        
        #Initialize Physics locations
        self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
        self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
        fluid['pore.'+self.name] = False
        fluid['pore.'+self.name][pores] = True
        fluid['throat.'+self.name] = False
        fluid['throat.'+self.name][throats] = True
        
    def physics_health(self):
        r'''
        Perform a check to find pores with overlapping or undefined Physics
        '''
        phys = self._fluid.physics()
        temp = sp.zeros((self._fluid.Np,))
        for item in phys:
            ind = self._fluid['pore.'+item]
            temp[ind] = temp[ind] + 1
        health = {}
        health['overlaps'] = sp.where(temp>1)[0].tolist()
        health['undefined'] = sp.where(temp==0)[0].tolist()
        return health
        
if __name__ == '__main__':
    print('none yet')


