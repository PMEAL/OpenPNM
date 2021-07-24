from openpnm.phases import GenericPhase
from openpnm.models.collections.phase import water


class Water(GenericPhase):
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
        self['pore.diffusivity'] = 1e-9
        self.models.update(water)
        self.regenerate_models()
