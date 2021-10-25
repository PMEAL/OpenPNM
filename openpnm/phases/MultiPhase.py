import numpy as np
import openpnm.models.misc as misc
from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class MultiPhase(GenericPhase):
    r"""
    Creates Phase object that represents a multiphase system consisting of
    a given list of OpenPNM Phase objects.

    Parameters
    ----------
    network : GenericNetwork
        The network to which this phase object will be attached.
    project : Project, optional
        The Project with which this phase should be associted. If a
        ``network`` is given then this is ignored and the Network's
        project is used. If a ``network`` is not given then this is
        mandatory.
    name : str, optional
        The name of the phase. This is useful to keep track of the objects
        throughout the simulation. The name must be unique to the project.
        If no name is given, one is generated.

    Examples
    --------
    >>> import scipy as sp
    >>> import openpnm as op
    >>> from openpnm.phases import Air, Water, MultiPhase

    >>> net = op.network.Cubic(shape=[5, 5, 5])
    >>> air = Air(network=net, name='air')  # Create two pure phases
    >>> water = Water(network=net, name='water')
    >>> mphase = MultiPhase(network=net, phases=[air, water], name='multi')
    >>> Ps = net['pore.coords'][:, 0] < 3  # Pick some pores to be air filled
    >>> Ts = net.find_neighbor_throats(pores=Ps)  # Find neighboring throats
    >>> Ts = net.tomask(throats=Ts)  # Convert throat indices to mask
    >>> mphase.set_occupancy(phase=air, Pvals=Ps, Tvals=Ts)  # Assign occupancies
    >>> mphase.set_occupancy(phase=water, Pvals=~Ps, Tvals=~Ts)

    Air and water have uniform viscosity values throughout, so the
    peak-to-peak distance is 0, while the mixture phase has the viscosity
    of air in some locations and water in others, hence a heterogenous
    viscosity array:

    >>> np.ptp(water['pore.viscosity']) == 0.
    True
    >>> np.ptp(air['pore.viscosity']) == 0.
    True
    >>> np.ptp(mphase['pore.viscosity']) == 0.
    False

    """

    def __init__(self, phases=[], settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(
            {
                'phases': [],
                'throat_occupancy': 'manual',
                'partition_coef_prefix': 'throat.partition_coef'
            }
        )
        self.settings.update(settings)

        self['pore.occupancy.all'] = np.zeros(self.Np, dtype=float)
        self['throat.occupancy.all'] = np.zeros(self.Nt, dtype=float)

        # Pressure/temperature must come from constituent phases
        self.pop('pore.temperature', None)
        self.pop('pore.pressure', None)

        # Add a dummy model for assembling the global parition
        # coefficients array. The dummy model only accepts 'target' and
        # kwargs. When a binary partition coefficient is added to the
        # MultiPhase, we manually append the binary partition coefficient
        # propname as a keyword argument to this dummy model, so that our
        # dependency handler is notified of this dependency!
        partition_coef_global = self.settings['partition_coef_prefix'] + ".all"
        self.add_model(propname=partition_coef_global, model=_dummy)

        # Add supplied phases to phases dict and initialize occupancy to 0
        self.add_phases(phases)

        logger.warning('MultiPhase is a beta feature. Functionality may change!')

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            vals = self.interleave_data(key)
        return vals

    def _get_phases(self):
        phases = {self.project[item].name: self.project[item]
                  for item in self.settings['phases']}
        return phases

    phases = property(fget=_get_phases)

    def _update_occupancy(self):
        r"""Updates 'occupancy.all' by summing up all occupancies"""
        for elem in ["pore", "throat"]:
            dict_ = self[f"{elem}.occupancy"]
            dict_.pop(f"{elem}.occupancy.all")
            self[f"{elem}.occupancy.all"] = np.sum(list(dict_.values()), axis=0)

    def add_phases(self, phases):
        r"""
        Adds supplied phases to the MultiPhase object and initializes
        occupancy to 0. It does not return a value.

        Parameters
        ----------
        phases : List[GenericPhase]

        """
        phases = np.array(phases, ndmin=1)
        for phase in phases:
            if phase.name not in self.settings['phases']:
                self.settings['phases'].append(phase.name)
                self[f'pore.occupancy.{phase.name}'] = 0.
                self[f'throat.occupancy.{phase.name}'] = 0.

    def set_binary_partition_coef(self, phases, model, **kwargs):
        r"""
        Sets binary partition coefficient as defined by the interface
        concentration ratio of phase 1 to phase 2.

        Parameters
        ----------
        phases : List[GenericPhase]
            List of the two phases for which the binary partition
            coefficient model is being added.
        model : OpenPNM model
            Model for calculating the binary partition coefficient.
        kwargs : dict
            Keyword arguments to be passed to the ``model``.

        Example
        -------
        >>> import openpnm as op
        >>> from openpnm.phases import Air, Water, MultiPhase
        >>> from openpnm.models.misc import constant

        >>> net = op.network.Cubic(shape=[5, 5, 5])
        >>> air = Air(network=net, name='air')  # Create two pure phases
        >>> water = Water(network=net, name='water')
        >>> mphase = MultiPhase(network=net, phases=[air, water], name='multi')

        Now, assign some pores to air and the rest to water

        >>> Ps = net['pore.coords'][:, 0] < 3  # Pick some pores to be air filled
        >>> Ts = net.find_neighbor_throats(pores=Ps)  # Find neighboring throats
        >>> Ts = net.tomask(throats=Ts)  # Convert throat indices to mask
        >>> mphase.set_occupancy(phase=air, Pvals=Ps, Tvals=Ts)  # Assign occupancies
        >>> mphase.set_occupancy(phase=water, Pvals=~Ps, Tvals=~Ts)

        Now, add an interface model for binary partition coefficient

        >>> mphase.set_binary_partition_coef(phases=[air, water],
        ...                                  model=constant,
        ...                                  value=0.5)

        Now, verify that K12 for interface throats is 0.5

        >>> Ts_interface = net.find_neighbor_throats(Ps, mode="xor")
        >>> K12_interface = mphase["throat.partition_coef.all"][Ts_interface]
        >>> assert K12_interface.mean() == 0.5

        Finally, verify that K12 for non-interface throats is 1.

        >>> Ts_rest = ~ np.isin(net.Ts, Ts_interface)
        >>> K12_rest = mphase["throat.partition_coef.all"][Ts_rest]
        >>> assert K12_rest.mean() == 1.

        """
        if np.size(phases) != 2:
            raise Exception("'phases' must contain exactly two elements!")

        # Add partition coefficient interface model to the MultiPhase
        propname_prefix = self.settings["partition_coef_prefix"]
        self._add_interface_prop(propname_prefix, phases, model, **kwargs)

        # Update global partition coef. model's args (for dependency handling)
        partition_coef_global = self.settings['partition_coef_prefix'] + ".all"
        K_global_model = self.models[partition_coef_global]
        # Append binary partition coefficient propname to model args
        propname = self._format_interface_prop(propname_prefix, phases)
        key = f"K_{phases[0].name}_{phases[1].name}"
        K_global_model[key] = propname

        # Regenerate 'throat.parition_coef.all'
        self.regenerate_models(partition_coef_global)

    def _add_interface_prop(self, propname, phases, model, **kwargs):
        r"""
        Helper method used to add interface models to the ``MultiPhase``
        object by augmenting.

        Notes
        -----
        For convenience and as an OpenPNM convention, the ``propname`` is
        augmented as outlined by the following example.

        Example
        -------
        Let's say the two phases corresponding to the interface model are
        named: 'air' and 'water', and the interface propname to be added
        is 'throat.foo'. After augmentation, 'throat.foo.air:water' will
        be the propname that's stored on the ``MultiPhase`` object.

        Note that because of this convention, the order of the phases that
        are passed to this method is important.

        """
        # Add "throat" keyword to the begining of propname if no identifier is found
        if propname.split(".")[0] not in ["pore", "throat"]:
            propname = f"throat.{propname}"
        # Check propname is throat property
        if propname.startswith("pore"):
            raise Exception("'propname' must be a throat property")
        # Add model to Multiphase
        propname = self._format_interface_prop(propname, phases)
        self.add_model(propname, model, **kwargs)

    def _format_interface_prop(self, propname, phases):
        r"""
        Returns formatted interface propname.

        Parameters
        ----------
        propname : str
            Dictionary key of the interface property.
        phases : List[GenericPhase]
            List of two phases that together form the interface property.

        Example
        -------
        >>> import openpnm as op
        >>> net = op.network.Cubic(shape=[5, 5, 5])
        >>> air = op.phases.Air(network=net, name='air')
        >>> water = op.phases.Water(network=net, name='water')
        >>> mphase = op.phases.MultiPhase(network=net, phases=[air, water])
        >>> mphase._format_interface_prop(propname="throat.foo", phases=[air, water])
        'throat.foo.air:water'

        """
        prefix = propname
        suffix = ":".join(phase.name for phase in phases)
        return f"{prefix}.{suffix}"

    def _get_phases_names(self, formatted_propname):
        r"""Retrieves phases' names from a formatted propname"""
        return formatted_propname.split(".")[2].split(":")

    def _assemble_partition_coef_global(self):
        r"""Updates the global partition coefficient array"""
        conns = self.network.conns
        partition_coef_models = []
        prefix = self.settings["partition_coef_prefix"]
        K_global = np.ones(self.Nt)
        partition_coef_global = prefix + ".all"

        # Find all binary partition coefficient models
        for model_key in self.models.keys():
            if model_key.startswith(prefix):
                if model_key != partition_coef_global:
                    partition_coef_models.append(model_key)

        # Modify the global partition coefficient for each phase pair
        for propname in partition_coef_models:
            K12 = self[propname]
            phases_names = self._get_phases_names(propname)
            S1 = self[f"pore.occupancy.{phases_names[0]}"][conns]
            S2 = self[f"pore.occupancy.{phases_names[1]}"][conns]
            mask = (S1[:, 0] + S2[:, 1]) == 2.
            K_global[mask] = K12[mask]
            mask = (S2[:, 0] + S1[:, 1]) == 2.
            K_global[mask] = 1 / K12[mask]

        return K_global

    def interleave_data(self, prop):
        r"""
        Gathers property values from component phases to build a single
        array.

        If the requested ``prop`` is not on this MultiPhase, then a search
        is conducted on all associated phase objects, and values from each
        are assembled into a single array.

        Parameters
        ----------
        prop : str
            The property to be retrieved.

        Returns
        -------
        ndarray
            An array containing the specified property retrieved from each
            component phase and assembled based on the specified mixing
            rule.

        """
        element = self._parse_element(prop)[0]
        vals = np.zeros([self._count(element=element)], dtype=float)
        # Retrieve property from constituent phases (weight = occupancy)
        try:
            for phase in self.phases.values():
                vals += phase[prop] * self[f"{element}.occupancy.{phase.name}"]
        # Otherwise - if not found - retrieve from super class
        except KeyError:
            vals = super().interleave_data(prop)

        # Check for consistency of occupancy values (i.e. add up to 1)
        if np.any(self[f"{element}.occupancy.all"] != 1.0):
            self._update_occupancy()
            if np.any(self[f"{element}.occupancy.all"] != 1.0):
                raise Exception(f"Occupancy doesn't add to unity in all {element}s")

        return vals

    def regenerate_models(self, propnames=None, exclude=[], deep=False):
        r"""
        Regenerate models associated with the Multiphase object

        This method works by first regenerating the models associated with
        the constituent phases, and then regenerating Multiphase models.

        Parameters
        ----------
        propnames : List[str] or str
            The list of property names to be regenerated. If None are
            given then ALL models are re-run (except for those whose
            ``regen_mode`` is 'constant').
        exclude : List[str]
            Since the default behavior is to run ALL models, this can be
            used to exclude specific models. It may be more convenient to
            supply as list of 2 models to exclude than to specify 8 models
            to include.
        deep : bool
            Specifies whether or not to regenerate models on all
            associated objects. For instance, if ``True``, then all
            Physics models will be regenerated when method is called on
            the corresponding Phase. The default is ``False``. The method
            does not work in reverse, so regenerating models on a Physics
            will not update a Phase.

        """
        # Regenerate models associated with phases within MultiPhase object
        for phase in self.phases.values():
            phase.regenerate_models(propnames=propnames, exclude=exclude, deep=deep)
        # Regenerate models specific to MultiPhase object
        super().regenerate_models(propnames=propnames, exclude=exclude, deep=deep)

    def _set_automatic_throat_occupancy(self, mode="mean"):
        r"""
        Automatically interpolates throat occupancy based on that in
        adjacent pores. This method doesn't return any value.

        Parameters
        ----------
        mode : str
            Interpolation method, ex. 'mean' sets the throat occupancy as
            the average of that in adjacent pores, while 'min' sets it to
            the minimum of the two. Options are 'mean', 'min', and 'max'.

        """
        self.settings['throat_occupancy'] = 'automatic'
        for phase in self.phases.values():
            self.add_model(
                propname=f"throat.occupancy.{phase.name}",
                model=misc.from_neighbor_pores,
                prop=f"pore.occupancy.{phase.name}",
                mode=mode
            )

    def set_occupancy(self, phase, Pvals=[], Tvals=[], pores=[], throats=[]):
        r"""
        Specifies occupancy of a phase in each pore and/or throat. This
        method doesn't return any value.

        Parameters
        ----------
        phase : GenericPhase
            The phase whose occupancy is being specified.
        Pvals : array_like
            The volume fraction of ``phase`` in each pore. This array must
            be Np-long, except when ``pores`` is also specified, where in
            that case they must be of equal length. If a scalar is
            received, it is applied to all the pores. If nothing is passed,
            ``Pvals=1.0`` will be assumed.
        Tvals : array_like
            The volume fraction of ``phase`` in each throat. This array
            must be Nt-long, except when ``throats`` is also specified,
            where in that case they must be of equal length. If a scalar
            is received it is applied to all the throats. If nothing is
            passed, ``Tvals=1.0`` will be assumed.
        pores : array_like
            The location of pores whose occupancy is to be set.
        throats : array_like
            The location of throats whose occupancy is to be set.

        """
        # TODO: pores/throats could also be masks

        Pvals = np.array(Pvals, ndmin=1)
        Tvals = np.array(Tvals, ndmin=1)
        pores = np.array(pores, ndmin=1)
        throats = np.array(throats, ndmin=1)

        # Check for size consistency of the arguments
        if Pvals.size and pores.size:
            if Pvals.size != 1 and pores.size != 1:
                if Pvals.size != pores.size:
                    raise Exception("Pvals and pores must be the same size.")
        if Tvals.size and throats.size:
            if Tvals.size != 1 and throats.size != 1:
                if Tvals.size != throats.size:
                    raise Exception("Tvals and throats must be the same size.")

        # Check if the passed phase is already part of MultiPhase object
        if phase not in self.project:
            raise Exception(f"{phase.name} doesn't belong to this project")
        # Add the passed phase to MultiPhase object if not found
        if phase.name not in self.settings['phases']:
            self.add_phases(phase)

        # Check for value consistency of the arguments
        if np.any(Pvals > 1.0) or np.any(Pvals < 0.0):
            logger.critical('Received Pvals contain volume fractions outside '
                            + 'the range of 0 to 1')
        if np.any(Tvals > 1.0) or np.any(Tvals < 0.0):
            logger.critical('Received Tvals contain volume fractions outside '
                            + 'the range of 0 to 1')

        if Pvals.size and not pores.size:
            pores = self.pores()
        if Tvals.size and not throats.size:
            throats = self.throats()

        if pores.size:
            Pvals = Pvals if Pvals.size else 1.0
            self['pore.occupancy.' + phase.name][pores] = Pvals
        if throats.size:
            Tvals = Tvals if Tvals.size else 1.0
            self['throat.occupancy.' + phase.name][throats] = Tvals

        if self.settings["throat_occupancy"] == "automatic":
            self.regenerate_models(propnames=f"throat.occupancy.{phase.name}")

        self._update_occupancy()


def _dummy(target, **kwargs):
    r"""Dummy model to store propnames as kwargs for dependency handling"""
    return target._assemble_partition_coef_global()
