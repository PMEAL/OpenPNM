import openpnm.models.physics as mods
from openpnm.utils import get_model_collection


def basic(regen_mode='deferred', domain=None):
    return get_model_collection(collection=_basic,
                                regen_mode=regen_mode,
                                domain=domain)


_basic = {
    'throat.hydraulic_conductance': {
        'model': mods.hydraulic_conductance.generic_hydraulic,
    },
    'throat.diffusive_conductance': {
        'model': mods.diffusive_conductance.generic_diffusive,
    },
    'throat.entry_pressure': {
        'model': mods.capillary_pressure.washburn,
    },
}
