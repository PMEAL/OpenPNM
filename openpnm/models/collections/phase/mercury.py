import openpnm.models as mods


mercury = {
    'throat.contact_angle': {
        'model': mods.misc.constant,
        'value': 140,
    },
    'pore.density': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 14280.9,
        'm': -2.47004,
    },
    'pore.molar_density': {
        'model': mods.phase.density.mass_to_molar,
    },
    'pore.surface_tension': {
        'model': mods.misc.linear,
        'prop': 'pore.temperature',
        'b': 0.56254,
        'm': -0.00028,
    },
    'pore.thermal_conductivity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [3.98691, 0.0170967, -0.0000063862],
    },
    'pore.viscosity': {
        'model': mods.misc.polynomial,
        'prop': 'pore.temperature',
        'a': [0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12],
    },
}

# Thermophysical Properties of Materials for Nuclear Engineering:
# IAEA, Vienna, 2008. ISBN 978-92-0-106508-7:
