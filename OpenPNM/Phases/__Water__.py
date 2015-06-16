# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm


class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

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
    >>> water = OpenPNM.Phases.Water(network=pn)
    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self._generate()

    def _generate(self):
        self['pore.molecular_weight'] = 0.01802               # kg/mol
        self['pore.critical_pressure'] = 2.2064E7             # Pa
        self['pore.critical_temperature'] = 647.1             # K
        self['pore.critical_volume'] = 0.003106               # kg/m3
        self['pore.contact_angle'] = 110.0                    # Degree
        self.models.add(propname='pore.density',
                        model=fm.density.water)               # kg/m3
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.standard)      # mol/m3
        self['pore.diffusivity'] = 1e-9                       # m2/s
        self.models.add(propname='pore.surface_tension',
                        model=fm.surface_tension.water)       # N/m
        self.models.add(propname='pore.thermal_conductivity',
                        model=fm.thermal_conductivity.water)  # W/m.K
        self.models.add(propname='pore.vapor_pressure',       # Pa
                        model=fm.vapor_pressure.antoine,
                        A=10.1965, B=1730.63, C=-39.720)
        self.models.add(propname='pore.viscosity',
                        model=fm.viscosity.water)             # kg/m.s
