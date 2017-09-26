# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm

#TODO ALL
class WaterVapor(GenericPhase):
    r"""
    Creates Phase object with preset values for WaterVapor

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which this phase object will be associated.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.
    
    http://webbook.nist.gov/chemistry/fluid/

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Phases.WaterVapor(network=pn)
    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self._generate()

    def _generate(self):
        mmh2= 0.0020159            # kg/mol of hydrogen
        # Atomic Diffusion Volumes (could be wrong)
        advh2=7.07
        self['pore.state'] = 0
        self['pore.atomic_diffusion_volume'] = 12.7
        self['pore.molecular_weight'] = 0.01802               # kg/mol
        self['pore.critical_pressure'] = 2.2064E7             # Pa
        self['pore.critical_temperature'] = 647.1             # K
        self['pore.critical_volume'] = 0.003106               # kg/m3
        self.models.add(propname='pore.density',
                        model=fm.density.ideal_gas)               # kg/m3
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.ideal_gas)      # mol/m3
                        
        self.models.add(propname='pore.diffusivity',
                        model=fm.diffusivity.fuller,
                        MA=mmh2, MB=self['pore.molecular_weight'],
                        vA=self['pore.atomic_diffusion_volume'], vB=advh2)
        
        self.models.add(propname='pore.thermal_conductivity',
                        model=fm.thermal_conductivity.water)  # W/m.K
                        
        
        self['pore.viscosity'] = 13.0e-6                      # Pa.s
        '''
        self.models.add(propname='pore.viscosity',
                        model=fm.viscosity.water)             # kg/m.s
        '''
