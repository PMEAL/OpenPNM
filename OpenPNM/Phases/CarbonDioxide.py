# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm

class CarbonDioxide(GenericPhase):
    r"""
    Creates Phase object with preset values for Carbon Dioxide

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which this phase object will be associated.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> water = OpenPNM.Phases.CO2(network=pn)
    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self._generate()

    def _generate(self):
        self['pore.molecular_weight'] = .04401               # kg/mol
        self['pore.critical_pressure'] = 7.3773E6            # Pa
        self['pore.critical_temperature'] = 304.1282         # K   
        self['pore.density'] = 1101        
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.standard)      # mol/m3
        self.models.add(propname='pore.surface_tension',
                        model=fm.surface_tension.eotvos,
                        k = 2.1E-7)                            # N/m

        print('Warning - CO2 phase infomation is not complete')
