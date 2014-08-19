"""
module Physics
===============================================================================

"""
import OpenPNM
import scipy as sp

class GenericPhysics(OpenPNM.Base.Core):
    r"""
    Generic class to generate Physics objects  

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Physics should be attached
        
    phase : OpenPNM Phase object 
        The Phase object to which this Physics applies
    
    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.
    
    name : str, optional
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.  
    
    """

    def __init__(self,network,phase,pores=[],throats=[],name=None,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        #Append objects self for internal access
        self._net = network
        self.name = name
        
        #Append self to other objects
        network._physics.append(self)
        phase._physics.append(self)
        self._phases.append(phase)
        
        #Initialize Physics locations
        self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
        self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
        phase['pore.'+self.name] = False
        phase['pore.'+self.name][pores] = True
        phase['throat.'+self.name] = False
        phase['throat.'+self.name][throats] = True
        
    def check_physics_health(self):
        r'''
        Perform a check to find pores with overlapping or undefined Physics
        '''
        phase = self._phases[0]
        phys = phase.physics()
        temp = sp.zeros((phase.Np,))
        for item in phys:
                ind = phase['pore.'+item]
                temp[ind] = temp[ind] + 1
        health = {}
        health['overlaps'] = sp.where(temp>1)[0].tolist()
        health['undefined'] = sp.where(temp==0)[0].tolist()
        return health
        
if __name__ == '__main__':
    print('none yet')


