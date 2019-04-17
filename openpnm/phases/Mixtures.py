from openpnm.phases import GenericPhase
import openpnm.models as mods


class SalineWater(GenericPhase):
    r"""
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.01802
        self['pore.critical_pressure'] = 2.2064E7
        self['pore.critical_temperature'] = 647.1
        self['pore.critical_volume'] = 0.003106
        self['pore.contact_angle'] = 110.0
        self['pore.electrical_conductivity'] = 1e-15

        self['pore.permittivity'] = 78.303
        self['pore.diffusivity.Na'] = 1.33e-09
        self['pore.diffusivity.Cl'] = 2.03e-09
        self['pore.valence.Na'] = 1
        self['pore.valence.Cl'] = -1

        self['throat.permittivity'] = 78.303
        self['throat.diffusivity.Na'] = 1.33e-09
        self['throat.diffusivity.Cl'] = 2.03e-09
        self['throat.valence.Na'] = 1
        self['throat.valence.Cl'] = -1

        self.add_model(propname='pore.density',
                       model=mods.phases.density.water)
        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.standard)
        self.add_model(propname='pore.surface_tension',
                       model=mods.phases.surface_tension.water)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.phases.thermal_conductivity.water)
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phases.vapor_pressure.antoine,
                       A=8.088, B=1750.71, C=236.191)
        self.add_model(propname='pore.viscosity',
                       model=mods.phases.viscosity.water)
        self.regenerate_models()
