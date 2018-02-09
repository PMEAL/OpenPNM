# -*- coding: utf-8 -*-
from openpnm.phases import GenericPhase
import openpnm.models.phase as fm


class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values for
    mercury.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    Notes
    -----
    The initial properties are all at std conditions of T = 298 K and P = 1
    atm.

    References
    ----------
    [1] Thermophysical Properties of Materials for Nuclear Engineering: IAEA,
        Vienna, 2008. ISBN 978-92-0-106508-7:

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> hg = op.phases.Mercury(network=pn)

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self['pore.molecular_weight'] = 0.2006                         # kg/mol
        self['pore.critical_pressure'] = 1.662E8                       # Pa
        self['pore.critical_temperature'] = 1733                       # K
        self['pore.critical_volume'] = 0.000189                        # kg/m3
        self['pore.contact_angle'] = 140                               # Degree
        self.add_model(propname='pore.vapor_pressure',                 # Pa
                       model=fm.vapor_pressure.antoine,
                       A=9.85767, B=3007.129, C=-10.001)
        self.add_model(propname='pore.density',                        # kg/m3
                       model=fm.misc.linear,
                       prop='pore.temperature',
                       b=14280.9, m=-2.47004)
        self.add_model(propname='pore.molar_density',
                       model=fm.molar_density.standard)                # mol/m3
        self.add_model(propname='pore.surface_tension',                # N/m
                       model=fm.misc.linear,
                       prop='pore.temperature',
                       b=0.56254, m=-0.00028)
        self.add_model(propname='pore.thermal_conductivity',           # W/m.K
                       model=fm.misc.polynomial,
                       prop='pore.temperature',
                       a=[3.98691, 0.0170967, -0.0000063862])
        self.add_model(propname='pore.viscosity',                      # kg/m.s
                       model=fm.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12])
        self.regenerate_models()
