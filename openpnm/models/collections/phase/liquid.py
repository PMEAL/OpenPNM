import openpnm.models as mods
from openpnm.utils import get_model_collection


def standard_liquid(regen_mode=None, domain=None):
    return get_model_collection(collection=_pure_liquid,
                                regen_mode=regen_mode,
                                domain=domain)


def standard_liquid_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_liquid_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_pure_liquid = {
    'pore.density': {
        'model': mods.phase.density.liquid_pure,
    },
}


_liquid_mixture = {
    'pore.density': {
        'model': mods.phase.density.liquid_mixture,
    },
}
