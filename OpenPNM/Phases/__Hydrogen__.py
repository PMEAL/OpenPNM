# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm


class Hydrogen(GenericPhase):
    r"""
    Creates Phase object with preset models and values for Hydrogen

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    Notes
    -----
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    References
    ----------
    

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> hydrogen = OpenPNM.Phases.Hydrogen(network=pn)

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self._generate()

    def _generate(self):
        # Molar Mass of water
        mmh20=0.01801528 # kg/mol
        # Atomic diffusion volumes (check these with Gostick?)
        advh20=12.7
        # Done
        self['pore.state'] = 0
        self['pore.atomic_diffusion_volume'] = 7.07
        self['pore.molecular_weight'] = 0.0020159            # kg/mol
        self['pore.critical_pressure'] = 1.30E6           # Pa
        self['pore.critical_temperature'] = 33.2          # K
        # WTF? TODO
        self['pore.critical_volume'] = 0.002917            # kg/m3
        # Useless
        self['pore.contact_angle'] = 110.0                 # Degree    
        self.models.add(propname='pore.density',
                        model=fm.density.ideal_gas)        # kg/m3
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.ideal_gas)  # mol/m3
        # Working on it
        self.models.add(propname='pore.diffusivity',
                        model=fm.diffusivity.fuller,
                        MA=mmh20, MB=self['pore.molecular_weight'],
                        vA=advh20, vB=self['pore.atomic_diffusion_volume'])
                        
        # I do not undertand this model
        '''
        self.models.add(propname='pore.thermal_conductivity',      # W/m.K
                        model=fm.misc.polynomial,
                        poreprop='pore.temperature',
                        # Where do these values come from?
                        a=[0.00422791, 0.0000789606, -1.56383E-08])
        '''
        
        # Thermal conductivity of hydrogen?
        self['pore.thermal_conductivity'] = 0.172                  # W/m.K
        # Viscosity of Hydrogen at 20C?
        self['pore.viscosity'] = 8.7454468E-6                      # kg/m.s
        # Wut
        '''
        self.models.add(propname='pore.viscosity',                 # kg/m.s
                        model=fm.misc.polynomial,
                        poreprop='pore.temperature',
                        a=[0.00000182082, 6.51815E-08, -3.48553E-11, 1.11409E-14])
        '''
