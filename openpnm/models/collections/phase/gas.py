import openpnm.models as mods
from openpnm.utils import get_model_collection


def standard_gas(regen_mode=None, domain=None):
    return get_model_collection(collection=_pure_gas,
                                regen_mode=regen_mode,
                                domain=domain)


def standard_gas_mixture(regen_mode=None, domain=None):
    return get_model_collection(collection=_gas_mixture,
                                regen_mode=regen_mode,
                                domain=domain)


_pure_gas = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
}


_gas_mixture = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
}
