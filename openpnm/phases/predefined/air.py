r"""

Example
-------
import openpnm as op
pn = op.network.Cubic(shape=[5, 1, 1])
p = op.phases.GenericPhase(network=pn)
p.models.update(op.phases.predefined.air)
p.regenerate_models()

"""

import openpnm.models as mods

air = {
       'pore.molecular_weight': {
           'model': mods.misc.constant,
           'value': 0.0291,
           'regen_mode': 'deferred',
           },
       'pore.critical_pressure': {
           'model': mods.misc.constant,
           'value': 3.786E6,
           'regen_mode': 'deferred',
           },
       'pore.critical_temperature': {
           'model': mods.misc.constant,
           'value': 132.5,
           'regen_mode': 'deferred',
           },
       'pore.critical_volume': {
           'model': mods.misc.constant,
           'value': 0.002917,
           'regen_mode': 'deferred',
           },
       'pore.contact_angle': {
           'model': mods.misc.constant,
           'value': 180.0,
           'regen_mode': 'deferred',
           },
       'pore.surface_tension': {
           'model': mods.misc.constant,
           'value': 0.072,
           'regen_mode': 'deferred',
           },
       'pore.molar_density': {
           'model': mods.phases.molar_density.ideal_gas,
           'regen_mode': 'deferred',
           },
       'pore.diffusivity': {
           'model': mods.phases.diffusivity.fuller,
           'MA': 0.032,
           'MB':0.028,
           'vA': 16.6,
           'vB': 17.9,
           'regen_mode': 'deferred',
           },
       'pore.thermal_conductivity': {
           'model': mods.misc.polynomial,
           'prop': 'pore.temperature',
           'a': [0.00422791, 0.0000789606, -1.56383E-08],
           'regen_mode': 'deferred',
           },
       'pore.electrical_conductivity': {
           'model': mods.misc.constant,
           'value': 1e-15,
           'regen_mode': 'deferred',
           },
       'pore.viscosity': {
           'model': mods.misc.polynomial,
           'prop': 'pore.temperature',
           'a': [0.00000182082, 6.51815E-08, -3.48553E-11,1.11409E-14],
           'regen_mode': 'deferred',
           },
       }
