from openpnm.phases import GenericPhase
import openpnm.models as mods


class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values and
    pore-scale models for mercury.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    project : OpenPNM Project object, optional
        Can be supplied instead of ``network``

    References
    ----------
    [1] Thermophysical Properties of Materials for Nuclear Engineering: IAEA,
        Vienna, 2008. ISBN 978-92-0-106508-7:

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> hg = op.phases.Mercury(network=pn)

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
    | 1 | pore.vapor_pressure  | model:           | antoine                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | A                | 9.85767                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | B                | 3007.129                 |
    +---+----------------------+------------------+--------------------------+
    |   |                      | C                | -10.001                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    | 2 | pore.density         | model:           | linear                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | b                | 14280.9                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | m                | -2.47004                 |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    | 3 | pore.molar_density   | model:           | standard                 |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | mol_weight       | pore.molecular_weight    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | density          | pore.density             |
    +---+----------------------+------------------+--------------------------+
    | 4 | pore.surface_tension | model:           | linear                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | b                | 0.56254                  |
    +---+----------------------+------------------+--------------------------+
    |   |                      | m                | -0.00028                 |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    | 5 | pore.thermal_cond... | model:           | polynomial               |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | a                | [3.98691, 0.0170967, ... |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    | 6 | pore.viscosity       | model:           | polynomial               |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | a                | [0.00355837, -1.00131... |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+

    """
    def __init__(self, name=None, **kwargs):
        super().__init__(name=name, **kwargs)

        self['pore.molecular_weight'] = 0.2006
        self['pore.critical_pressure'] = 1.662E8
        self['pore.critical_temperature'] = 1733
        self['pore.critical_volume'] = 0.000189
        self['pore.contact_angle'] = 140.0
        self['pore.electrical_conductivity'] = 1e6
        self['pore.diffusivity'] = 1e-15
        self.add_model(propname='pore.vapor_pressure',
                       model=mods.phases.vapor_pressure.antoine,
                       A=9.85767, B=3007.129, C=-10.001)
        self.add_model(propname='pore.density',
                       model=mods.misc.linear,
                       prop='pore.temperature',
                       b=14280.9, m=-2.47004)
        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.standard)
        self.add_model(propname='pore.surface_tension',
                       model=mods.misc.linear,
                       prop='pore.temperature',
                       b=0.56254, m=-0.00028)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[3.98691, 0.0170967, -0.0000063862])
        self.add_model(propname='pore.viscosity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00355837, -0.0000100131, 1.23684E-08, -5.1684E-12])
        self.regenerate_models()
