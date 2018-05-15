from openpnm.phases import GenericPhase
import openpnm.models as mods


class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which this phase object will be associated.

    Notes
    -----
    The initial properties are all at std conditions of T = 298 K and
    P = 1 atm.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> water = op.phases.Water(network=pn)
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.01802
        self['pore.critical_pressure'] = 2.2064E7
        self['pore.critical_temperature'] = 647.1
        self['pore.critical_volume'] = 0.003106
        self['pore.contact_angle'] = 110.0
        self['pore.diffusivity'] = 1e-9
        self.add_model(propname='pore.density',
                       model=mods.phases.density.water)
        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.standard)
        self.add_model(propname='pore.surface_tension',
                       model=mods.phases.surface_tension.water)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.phases.thermal_conductivity.water)
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phases.vapor_pressure.antoine,
                       A=8.088, B=1750.71, C=236.191)
        self.add_model(propname='pore.viscosity',
                       model=mods.phases.viscosity.water)
        self.regenerate_models()
