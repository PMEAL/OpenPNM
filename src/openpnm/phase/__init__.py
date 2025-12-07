r"""
Phase
=====

This module contains classes generating and handling thermophysical properties

"""


def _fetch_chemical_props(a):
    temp = {}
    k = 1.380649e-23  # Boltzmann constant
    if not hasattr(a, 'zs'):  # If NOT a mixture
        temp['CAS'] = a.CAS
        temp['common_name'] = a.name
        temp['charge'] = a.charge
        temp['formula'] = a.formula
        temp['boiling_temperature'] = a.Tb
        temp['melting_temperature'] = a.Tm
        temp['triple_point_temperature'] = a.Tt
        temp['triple_point_pressure'] = a.Pt
        temp['dipole_moment'] = a.dipole
        temp['LJ_diameter'] = a.molecular_diameter
        temp['LJ_energy'] = a.Stockmayer*k
        # Change temperature to Tb to evaluate surface tension at Tb
        a.T = a.Tb
        temp['surface_tension_Tb'] = a.sigma
        temp['molar_volume_Tb'] = a.Vml_Tb  # Already at Tb
    else:
        temp['formula'] = a.formulas
        temp['mole_fractions'] = a.zs
        temp['CAS'] = None
    temp['molecular_weight'] = a.MW
    temp['critical_temperature'] = a.Tc
    temp['critical_pressure'] = a.Pc
    temp['critical_volume'] = a.Vc
    temp['critical_compressibilty_factor'] = a.Zc
    temp['acentric_factor'] = a.omega
    return temp


def _get_mixture_model_args(
    phase,
    composition='xs',
    args={
        'mus': 'pore.viscosity',
        'MWs': 'param.molecular_weight',
    }
):
    from openpnm.models.phase.misc import mole_to_mass_fraction
    from numpy import vstack
    vals = {}
    if composition in ['ws']:
        temp = vstack(list(mole_to_mass_fraction(phase=phase).values()))[:, 0]
        vals[composition] = temp
    else:
        temp = vstack(list(phase['pore.mole_fraction'].values()))[:, 0]
        vals[composition] = temp
    for item in args.keys():
        temp = vstack(list(phase.get_comp_vals(args[item]).values()))[:, 0]
        vals[item] = temp
    return vals


from ._phase import *
from ._mixture import *
from ._species import *
from ._air import *
from ._water import *
from ._mercury import *
