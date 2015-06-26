# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm


class Air(GenericPhase):
    r"""
    Creates Phase object with preset models and values for air

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    Notes
    -----
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    References
    ----------
    [1] E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity
    Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of Thermophysics,
    Vol. 25, No. 1, January 2004, pp. 21-69

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> air = OpenPNM.Phases.Air(network=pn)

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self._generate()

    def _generate(self):
        self['pore.molecular_weight'] = 0.02896            # kg/mol
        self['pore.critical_pressure'] = 3.786E6           # Pa
        self['pore.critical_temperature'] = 132.5          # K
        self['pore.critical_volume'] = 0.002917            # kg/m3
        self['pore.contact_angle'] = 110.0                 # Degree
        self.models.add(propname='pore.density',
                        model=fm.density.ideal_gas)        # kg/m3
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.ideal_gas)  # mol/m3
        self.models.add(propname='pore.diffusivity',
                        model=fm.diffusivity.fuller,
                        MA=0.032, MB=0.028,
                        vA=16.6, vB=17.9)
        self.models.add(propname='pore.thermal_conductivity',      # W/m.K
                        model=fm.misc.polynomial,
                        poreprop='pore.temperature',
                        a=[0.00422791, 0.0000789606, -1.56383E-08])
        self.models.add(propname='pore.viscosity',                 # kg/m.s
                        model=fm.misc.polynomial,
                        poreprop='pore.temperature',
                        a=[0.00000182082, 6.51815E-08, -3.48553E-11, 1.11409E-14])
