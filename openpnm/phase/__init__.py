def _fetch_chemical_props(a):
    temp = {}
    if not hasattr(a, 'zs'):
        temp['CAS'] = a.CAS
        temp['common_name'] = a.name
        temp['charge'] = a.charge
        temp['formula'] = a.formula
        temp['boiling_temperature'] = a.Tb
        temp['melting_temperature'] = a.Tm
        temp['triple_point_temperature'] = a.Tt
        temp['triple_point_pressure'] = a.Pt
        temp['dipole_moment'] = a.dipole
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


from ._phase import Phase
from ._air import Air
from ._water import Water
from ._mercury import Mercury
from ._mixture import *
from ._species import *
