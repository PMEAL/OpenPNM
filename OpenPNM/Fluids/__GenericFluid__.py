# -*- coding: utf-8 -*-

import OpenPNM
import scipy as sp

class GenericFluid(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericFluid - Base class to manipulate fluid properties

    Parameters
    ----------


    Examples
    --------
    >>> print 'nothing yet'

    .. note::
    n/a

    """    
    def __init__(self,fluid_dict,**logvars):
        r"""
        Initialize
        """
        super(GenericFluid,self).__init__(**logvars)
        self.indent = ""
        self._logger.debug("Construct class")
        self._fluid_dict = fluid_dict
        self._property_list = ['name','diffusivity','viscosity','molardensity']
        
    def update(self,pn=OpenPNM.Network.GenericNetwork()):
        r"""
        step through all values in the property_list and create pore and throat conditions accordingly
        """
        phase = self._find_phase(pn)
        self._update_name(pn,phase)
        self._update_diffusivity(pn,phase)
        self._update_viscosity(pn,phase)
        self._update_molardensity(pn,phase)
        
    def add_property(self,prop_name,value):
        if prop_name in self._property_list:
            self._fluid_dict[prop_name] = value
        else:
            self._logger.error("'"+prop_name+"'"+' is not a valid fluid property! \n\t\tValid names are:'+str(self._property_list.__str__()))
        
    def _find_phase(self,pn=OpenPNM.Network.GenericNetwork()):
        r"""
        determine whether the phase for this fluid has been specified
        if 'phase' is within the fluid_dict, then use that number, overwriting anything that is already there
        if 'phase' is absent, check the pore network for that phase, if it isn't there, create a new phase, starting with '1'
        """
        if 'phase' in self._fluid_dict.keys():
            phase = self._fluid_dict['phase']
        else:
            if self._fluid_dict['name'] in pn.phases:
                for i in range(length(pn.phases)):
                    if pn.phases[i] == self._fluid_dict['name']:
                        phase = i+1
            else:
                phase = sp.size(pn.phases)+1
        return phase
    
    def _update_name(self,pn,phase):
        r"""
        Updates the fluid name

        Notes
        -----
        This method is not implemented in the GenericFluid, and must be subclassed in each algorithm as necessary
        """
        self._logger.error("_update_name: not implemented")
    def _update_diffusivity(self,pn,phase):
        r"""
        Updates the fluid diffusivity

        Notes
        -----
        This method is not implemented in the GenericFluid, and must be subclassed in each algorithm as necessary
        """
        self._logger.error("_update_diffusivity: not implemented")
    def _update_viscosity(self,pn,phase):
        r"""
        Updates the fluid viscosity

        Notes
        -----
        This method is not implemented in the GenericFluid, and must be subclassed in each algorithm as necessary
        """
        self._logger.error("_update_viscosity: not implemented")
    def _update_molardensity(self,pn,phase):
        r"""
        Updates the fluid molardensity

        Notes
        -----
        This method is not implemented in the GenericFluid, and must be subclassed in each algorithm as necessary
        """
        self._logger.error("_update_molardensity: not implemented")
        
if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)