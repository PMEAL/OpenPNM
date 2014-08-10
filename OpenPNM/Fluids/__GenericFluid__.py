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

class GenericFluid(OpenPNM.Core):
    r"""
    Base class to generate a generic fluid object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common fluids are included with OpenPNM and can be found under OpenPNM.Fluids.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this fluid should be attached
    name : str
        A unique string name to identify the fluid object, typically same as 
        instance name but can be anything.
    init_cond : dictionary, optional
        A dictionary of 'key':value pairs to initize fluid properties.  If omitted
        the temperature and pressure of the fluid are set to STP.  A typical 
        additiona would be 'mole_fraction' for a mixture.  Temperature and
        pressure can be set through init_cond, but will be set to STP if not 
        only other properties are included.
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """
    def __init__(self,network,name=None,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        # Attach objects to self for internal access
        self._net = network
        self.name = name
        
        # Append this Fluid to the Network
        network._fluids.append(self) 
        
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
        super(GenericFluid,self).__setitem__(prop,value)
        
    def __getitem__(self,key):
        if key not in self.keys():
            self._logger.debug(key+' not on Fluid, constructing data from Physics')
            return self.interleave_data(key,sources=self._physics)
        else:
            return super().__getitem__(key)
    
if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    fluid = OpenPNM.Fluids.GenericFluid(name='test_fluid',network=pn)
    import doctest
    doctest.testmod(verbose=True)



