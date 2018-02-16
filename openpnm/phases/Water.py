from openpnm.phases import GenericPhase
from openpnm.phases import models as fm


class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which this phase object will be associated.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> water = op.phases.Water(network=pn)
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.01802               # kg/mol
        self['pore.critical_pressure'] = 2.2064E7             # Pa
        self['pore.critical_temperature'] = 647.1             # K
        self['pore.critical_volume'] = 0.003106               # kg/m3
        self['pore.contact_angle'] = 110.0                    # Degree
        self['pore.diffusivity'] = 1e-9                       # m2/s
        self.add_model(propname='pore.density',
                       model=fm.density.water)               # kg/m3
        self.add_model(propname='pore.molar_density',
                       model=fm.molar_density.standard)      # mol/m3
        self.add_model(propname='pore.surface_tension',
                       model=fm.surface_tension.water)       # N/m
        self.add_model(propname='pore.thermal_conductivity',
                       model=fm.thermal_conductivity.water)  # W/m.K
        self.add_model(propname='pore.vapor_pressure',       # Pa
                       model=fm.vapor_pressure.antoine,
                       A=8.088, B=1750.71, C=236.191)
        self.add_model(propname='pore.viscosity',
                       model=fm.viscosity.water)             # kg/m.s
        self.regenerate_models()
