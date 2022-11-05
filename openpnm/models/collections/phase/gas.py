import openpnm.models as mods


standard_gas = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity_gas': {
        'model': mods.phase.heat_capacity.gas_pure_TRC,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_pure_TRC,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_pure_gismr,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_pure_gesmr,
    },
}


standard_gas_mixture = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_mixture_yweighted,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_mixture_whz,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_mixture_hz,
    },
}

binary_gas_mixture = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.heat_capacity': {
        'model': mods.phase.heat_capacity.gas_mixture_yweighted,
    },
    'pore.thermal_conductivity': {
        'model': mods.phase.thermal_conductivity.gas_mixture_whz,
    },
    'pore.viscosity': {
        'model': mods.phase.viscosity.gas_mixture_hz,
    },
    'pore.diffusivity': {
        'model': mods.phase.diffusivity.gas_mixture_ce,
    },
}
