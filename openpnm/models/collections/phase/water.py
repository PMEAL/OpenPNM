import openpnm.models as mods

water = {
    'pore.density': {
        'model': mods.phases.density.water,
        },
    'pore.molar_density': {
        'model': mods.phases.molar_density.standard,
        },
    'pore.surface_tension': {
        'model': mods.phases.surface_tension.water,
        },
    'pore.thermal_conductivity': {
        'model': mods.phases.thermal_conductivity.water,
        },
    'pore.vapor_pressure': {
        'model': mods.phases.vapor_pressure.antoine,
        'A': 8.088,
        'B': 1750.71,
        'C': 236.191,
        },
    'pore.viscosity': {
        'model': mods.phases.viscosity.water,
        },
    }
