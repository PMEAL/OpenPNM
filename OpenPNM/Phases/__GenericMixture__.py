"""
module Physics
===============================================================================

"""
import sys, os
import OpenPNM
import scipy as sp

class GenericMixture(OpenPNM.Core):
    r"""
    Base class to generate a generic phase object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common phases are included with OpenPNM and can be found under OpenPNM.Phases.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this phase should be attached
    components : list of OpenPNM Phase objects
        The phase objects which the mixture consists of.
    name : str, optional
        A unique string name to identify the phase object, typically same as 
        instance name but can be anything.

    """
    def __init__(self,network,components=[],name=None,**kwargs):
        super(GenericMixture,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        # Attach objects to self for internal access
        self._net = network
        self.name = name
        [self._phases.append(comp) for comp in components]
        
        # Append this Phase to the Network
        network._phases.append(self)
        
        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']        
        
        # Set default T and P since most propery models require it
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        
    def __setitem__(self,prop,value):
        for phys in self._physics:
            if prop in phys.keys():
                self._logger.error(prop+' is already defined in at least one associated Geometry object')
                return
        super(GenericMixture,self).__setitem__(prop,value)
        
    def __getitem__(self,key):
        if key not in self.keys():
            self._logger.debug(key+' not on Phase, constructing data from Physics')
            return self._interleave_data(key,sources=self._physics)
        else:
            return super(GenericMixture,self).__getitem__(key)
            
    def set_component(self,phase,mode='add'):
        r'''
        '''
        if type(phase) == str:
            phase = self.find_object(obj_name=phase)
        if mode == 'add':
            self._phases.append(phase)
        elif mode == 'remove':
            self._phases.remove(phase)
            
    def check_mixture_health(self):
        r'''
        '''
        mole_sum = sp.zeros((self.Np,))
        for phase in self._phases:
            try:
                mole_sum = mole_sum + phase['pore.mole_fraction']
            except:
                pass
        return mole_sum
    
if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    phase = OpenPNM.Phases.GenericPhase(name='test_phase',network=pn)
    import doctest
    doctest.testmod(verbose=True)

