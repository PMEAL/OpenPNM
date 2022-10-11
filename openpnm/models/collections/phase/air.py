import openpnm.models as mods


air = {
    'pore.density': {
        'model': mods.phase.density.ideal_gas,
    },
    'pore.molar_density': {
        'model': mods.phase.density.mass_to_molar,
    },
    'pore.diffusivity': {
        'model': mods.phase.diffusivity.gas_mixture_fesg,
        'MWs': [31.9988, 28.0134],
        'Vdms': [16.6, 17.9],
    },
    'pore.thermal_conductivity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00422791, 0.0000789606, -1.56383E-08],
    },
    'pore.viscosity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00000182082, 6.51815E-08, -3.48553E-11, 1.11409E-14],
    },
}
