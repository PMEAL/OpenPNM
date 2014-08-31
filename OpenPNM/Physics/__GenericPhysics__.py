"""
===============================================================================
module __Physics__: Base class for mananging pore-scale Physics properties
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

    def __init__(self,network=None,phase=None,pores=[],throats=[],name=None,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")

        #Initialize locations
        self['pore.all'] = sp.array([],dtype=bool)
        self['throat.all'] = sp.array([],dtype=bool)

        #Associate with Network
        if network == None:
            self._net = OpenPNM.Network.GenericNetwork()
        else:
            self._net = network  # Attach network to self
            self._net._physics.append(self)  # Register self with network
        self.name = name

        #Associate with Phase
        if phase == None:
            self._phases.append(OpenPNM.Phases.GenericPhase())
        else:
            phase._physics.append(self)  # Register self with phase
            self._phases.append(phase)  # Register phase with self

        #Initialize a label dictionary in the associated fluid
        self._phases[0]['pore.'+self.name] = False
        self._phases[0]['throat.'+self.name] = False
        self.set_locations(pores=pores,throats=throats)

    def set_locations(self,pores=[],throats=[]):
        r'''
        This method can be used to set the pore and throats locations of an
        *empty* object.  Once locations have been set they can not be changed.
        
        Parameters
        ----------
        pores and throats : array_like
            The list of pores and/or throats where the object should be applied.
            
        Notes
        -----
        This method is intended to assist in the process of loading saved
        objects.  Save data can be loaded onto an empty object, then the object 
        can be reassociated with a Network manually by setting the pore and 
        throat locations on the object.  
        '''
        if pores != []:
            #Initialize locations
            self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
            self['pore.map'] = pores
            #Specify Physics locations in Phase dictionary
            self._phases[0]['pore.'+self.name][pores] = True
        if throats != []:
            #Initialize locations
            self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
            self['throat.map'] = throats
            #Specify Physics locations in Phase dictionary
            self._phases[0]['throat.'+self.name][throats] = True

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


