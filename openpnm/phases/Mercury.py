from openpnm.phases import GenericPhase
from openpnm.models.collections.phase import mercury


class Mercury(GenericPhase):
    r"""


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
        self.models.update(mercury)
        self.regenerate_models()
