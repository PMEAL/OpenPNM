from openpnm.phases import GenericPhase
import openpnm.models as mods


class Air(GenericPhase):
    r"""
    Creates Phase object with preset models and values for air.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project is
        used.  If a ``network`` is not given then this is mandatory.

    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.  If
        no name is given, one is generated.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> air = op.phases.Air(network=pn)

    Notes
    -----
    The table below shows all of the pore-scale models that are included with
    this class to calculate the physical properties of this fluid as functions
    of the relevant state variables.

    This object is initialized at standard conditions of 298 K and 101325 Pa.
    If these conditions are changed, the dependent properties can be
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
    | 1 | pore.molar_density   | model:           | ideal_gas                |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | pressure         | pore.pressure            |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    | 2 | pore.diffusivity     | model:           | fuller                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | MA               | 0.032                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | MB               | 0.028                    |
    +---+----------------------+------------------+--------------------------+
    |   |                      | vA               | 16.6                     |
    +---+----------------------+------------------+--------------------------+
    |   |                      | vB               | 17.9                     |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    |   |                      | temperature      | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | pressure         | pore.pressure            |
    +---+----------------------+------------------+--------------------------+
    | 3 | pore.thermal_cond... | model:           | polynomial               |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | a                | [0.00422791, 7.89606e... |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+
    | 4 | pore.viscosity       | model:           | polynomial               |
    +---+----------------------+------------------+--------------------------+
    |   |                      | prop             | pore.temperature         |
    +---+----------------------+------------------+--------------------------+
    |   |                      | a                | [1.82082e-06, 6.51815... |
    +---+----------------------+------------------+--------------------------+
    |   |                      | regen_mode       | normal                   |
    +---+----------------------+------------------+--------------------------+

    References
    ----------
    The pore scale models for this class are taken from:

    [1] E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity
    Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of Thermophysics,
    Vol. 25, No. 1, January 2004, pp. 21-69

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.0291
        self['pore.critical_pressure'] = 3.786E6
        self['pore.critical_temperature'] = 132.5
        self['pore.critical_volume'] = 0.002917
        self['pore.contact_angle'] = 180.0
        self['pore.surface_tension'] = 0.072
        self.add_model(propname='pore.molar_density',
                       model=mods.phases.molar_density.ideal_gas)
        self.add_model(propname='pore.diffusivity',
                       model=mods.phases.diffusivity.fuller,
                       MA=0.032, MB=0.028,
                       vA=16.6, vB=17.9)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00422791, 0.0000789606, -1.56383E-08])
        self.add_model(propname='pore.electrical_conductivity',
                       model=mods.misc.constant,
                       value=1e-15)
        self.add_model(propname='pore.viscosity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00000182082, 6.51815E-08, -3.48553E-11,
                          1.11409E-14])
