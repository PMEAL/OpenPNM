from openpnm.phases import GenericPhase
import openpnm.models as mods


class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values and
    pore-scale models for mercury.

    Parameters
    ----------
    %(GenericPhase.parameters)s

    References
    ----------
    The correlations and constants for this class were taken from:

    ::

        Thermophysical Properties of Materials for Nuclear Engineering:
        IAEA, Vienna, 2008. ISBN 978-92-0-106508-7:

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)

        self['pore.molecular_weight'] = 0.2006
        self['pore.critical_pressure'] = 1.662E8
        self['pore.critical_temperature'] = 1733
        self['pore.critical_volume'] = 0.000189
        self['pore.contact_angle'] = 140.0
        self['pore.electrical_conductivity'] = 1e6
        self['pore.diffusivity'] = 1e-15
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phases.vapor_pressure.antoine,
                       A=9.85767, B=3007.129, C=-10.001)
        self.add_model(propname='pore.density',
                       model=mods.misc.linear,
                       prop='pore.temperature',
                       b=14280.9, m=-2.47004)
        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.standard)
        self.add_model(propname='pore.surface_tension',
                       model=mods.misc.linear,
                       prop='pore.temperature',
                       b=0.56254, m=-0.00028)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[3.98691, 0.0170967, -0.0000063862])
        self.add_model(propname='pore.viscosity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12])
        self.regenerate_models()
