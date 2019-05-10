from openpnm.phases import GenericPhase
import openpnm.models as mods


class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which this phase object will be associated.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> water = op.phases.Water(network=pn)

    Notes
    -----
    The table below shows all of the pore-scale models that are included with
    this class to calculate the physical properties of this fluid as functions
    of the relevant state variables.

    This object is initialized at standard conditions of 298 K and 101325 Pa.
    If these conditions are changed the dependent properties can be
    recalculated by calling ``regenerate_models``.

    All of these parameters can be adjusted manually by editing the entries in
    the **ModelsDict** stored in the ``models`` attribute of the object.

    For a full listing of models and their parameters use ``print(obj.models)``
    where ``obj`` is the handle to the object.

    In addition to these models, this class also has a number of constant
    values assigned to it which can be found by running
    ``props(mode='constants')``.

    +---+----------------------+------------------+--------------------------+
    | # | Property Name        | Parameter        | Value                    |
    +===+======================+==================+==========================+
    | 1 | pore.density         | model:           | water                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | salinity         | pore.salinity            |
    +---+----------------------+------------------+--------------------------+
    | 2 | pore.molar_density   | model:           | standard                 |
    +---+----------------------+------------------+--------------------------+
    |   |                      | mol_weight       | pore.molecular_weight    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | density          | pore.density             |
    +---+----------------------+------------------+--------------------------+
    | 3 | pore.surface_tension | model:           | water                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | salinity         | pore.salinity            |
    +---+----------------------+------------------+--------------------------+
    | 4 | pore.thermal_cond... | model:           | water                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | salinity         | pore.salinity            |
    +---+----------------------+------------------+--------------------------+
    | 5 | pore.vapor_pressure  | model:           | antoine                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | A                | 8.088                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | B                | 1750.71                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | C                | 236.191                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    | 6 | pore.viscosity       | model:           | water                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | salinity         | pore.salinity            |
    +---+----------------------+------------------+--------------------------+


    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.01802
        self['pore.critical_pressure'] = 2.2064E7
        self['pore.critical_temperature'] = 647.1
        self['pore.critical_volume'] = 0.003106
        self['pore.contact_angle'] = 110.0
        self['pore.electrical_conductivity'] = 1e-15
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
