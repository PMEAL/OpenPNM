from openpnm.phase import GenericPhase
import openpnm.models as mods
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    %(GenericPhase.parameters)s

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.01802
        self['pore.critical_pressure'] = 2.2064E7
        self['pore.critical_temperature'] = 647.1
        self['pore.critical_volume'] = 0.003106
        self['pore.contact_angle'] = 110.0
        self['pore.electrical_conductivity'] = 1e-15
        self['pore.diffusivity'] = 1e-9

        self.add_model(propname='pore.density',
                       model=mods.phase.density.water)
        self.add_model(propname='pore.molar_density',
                       model=mods.phase.molar_density.standard)
        self.add_model(propname='pore.surface_tension',
                       model=mods.phase.surface_tension.water)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.phase.thermal_conductivity.water)
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phase.vapor_pressure.antoine,
                       A=8.088, B=1750.71, C=236.191)
        self.add_model(propname='pore.viscosity',
                       model=mods.phase.viscosity.water)
        self.regenerate_models()
