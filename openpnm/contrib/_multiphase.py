import logging

import numpy as np

import openpnm.models.misc as misc
from openpnm.phase import Phase as Phase
from openpnm.utils import Docorator, TypedSet

logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = [
    'MultiPhase',
    'multiphase_diffusion',
]


@docstr.dedent
class MultiPhaseSettings:
    """
    Parameters
    ----------
    %(PhaseSettings.parameters)s
    phases : list of strings
        The name of the phase objects which comprize this multiphase object
    throat_occupancy : str
        Indicates how throat occupancy is determined.  Options are:

        =========== ==========================================================
        mode        description
        =========== ==========================================================
        'automatic' The occupancy of each throat is calculated based on the
                    occupancy of the two neighboring pores
        'manual'    The occupancy of each throat must be set by hand by
                    assigning values to 'throat.occupancy.{phase.name}'
        =========== ==========================================================

    partition_coef_prefix : str
        The throat property which contains the partition coefficient values
    """
    phases = TypedSet(types=[str])
    throat_occupancy = 'manual'
    partition_coef_prefix = 'throat.partition_coef'


@docstr.dedent
class MultiPhase(Phase):
    """
    Creates a Phase object that represents a multiphase system consisting of
    a given list of Phases.

    Parameters
    ----------
    phases : list[Phase]
        A list containing the phase objects that make up the multiphase system
    %(Phase.parameters)s

    Notes
    -----
    This class assumes that only a SINGLE phase exists in each pore/throat.

    """

    def __init__(self, phases=[], name='mphase_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(MultiPhaseSettings())

        # Pressure/temperature must come from constituent phases
        self.pop('pore.temperature', None)
        self.pop('pore.pressure', None)

        # Initialize the partition coefficient, K
        self._K = np.ones(self.Nt, dtype=float)
        prefix = self.settings["partition_coef_prefix"]
        self[f"{prefix}.global"] = self._K

        # Add supplied phases to phases dict and initialize occupancy to 0
        self.add_phases(phases)

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            vals = self._interleave_data(key)
        return vals

    @property
    def phases(self):
        phases = {self.project[item].name: self.project[item]
                  for item in self.settings['phases']}
        return phases

    @property
    def K(self):
        self._build_K()
        return self._K

    @K.setter
    def K(self, value):
        self._K = value

    def add_phases(self, phases):
        """
        Adds supplied phases to MultiPhase object and sets occupancy to 0.

        Parameters
        ----------
        phases : list[Phase] or Phase

        """
        phases = np.array(phases, ndmin=1)
        for phase in phases:
            if phase.name in self.settings["phases"]:
                continue
            self.settings['phases'].add(phase.name)
            self[f'pore.occupancy.{phase.name}'] = 0.0
            self[f'throat.occupancy.{phase.name}'] = 0.0

    def set_occupancy(self, phase, *, pores=[], throats=[], values=1):
        r"""
        Specifies occupancy of a phase in each pore or throat. This
        method doesn't return any value.

        Parameters
        ----------
        phase : Phase
            The phase whose occupancy is being specified.
        pores : ndarray
            The location of pores whose occupancy is to be set.
        throats : ndarray
            The location of throats whose occupancy is to be set.
        values : ndarray or float
            Pore/throat occupancy values.

        """
        pores = np.array(pores, ndmin=1)
        throats = np.array(throats, ndmin=1)

        if not(pores.size ^ throats.size):
            raise Exception("Must either pass 'pores' or 'throats'")
        if phase not in self.project:
            raise Exception(f"{phase.name} doesn't belong to this project")
        self.add_phases(phase)

        if pores.size:
            self[f'pore.occupancy.{phase.name}'][pores] = values
        if throats.size:
            self[f'throat.occupancy.{phase.name}'][throats] = values

        if self.settings["throat_occupancy"] == "automatic":
            self.regenerate_models(propnames=f"throat.occupancy.{phase.name}")

    def regenerate_models(self, propnames=None, exclude=[]):
        r"""
        Regenerate models associated with the Multiphase object

        This method works by first regenerating the models associated with
        the constituent phases, and then regenerating Multiphase models.

        Parameters
        ----------
        propnames : list[str] or str
            The list of property names to be regenerated. If None are
            given then ALL models are re-run (except for those whose
            ``regen_mode`` is 'constant').
        exclude : list[str]
            Since the default behavior is to run ALL models, this can be
            used to exclude specific models. It may be more convenient to
            supply as list of 2 models to exclude than to specify 8 models
            to include.

        """
        # Regenerate models associated with phases within MultiPhase object
        for phase in self.phases.values():
            phase.regenerate_models(propnames=propnames, exclude=exclude)
        # Regenerate models specific to MultiPhase object
        super().regenerate_models(propnames=propnames, exclude=exclude)

    def set_binary_partition_coef(self, phases, model, **kwargs):
        """
        Sets binary partition coefficient as defined by the interface
        concentration ratio of phase 1 to phase 2.

        Parameters
        ----------
        phases : list[Phase]
            List of the two phases for which the binary partition
            coefficient model is being added.
        model : OpenPNM model
            Model for calculating the binary partition coefficient.
        kwargs : dict
            Keyword arguments to be passed to the ``model``.

        """
        assert len(phases) == 2
        # Add partition coefficient interface model to the MultiPhase
        propname_prefix = self.settings["partition_coef_prefix"]
        self._add_interface_prop(propname_prefix, phases, model, **kwargs)
        self._build_K()

    def _add_interface_prop(self, propname, phases, model, **kwargs):
        """
        Adds an interface model to the MultiPhase object. See Notes.

        Notes
        -----
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
        """Formats propname as {propname}.{phase[0].name}:{phase[1].name}"""
        prefix = propname
        suffix = ":".join(phase.name for phase in phases)
        return f"{prefix}.{suffix}"

    def _get_phase_labels(self, formatted_propname):
        """Retrieves phases names from a formatted propname"""
        assert ":" in formatted_propname
        return formatted_propname.split(".")[-1].split(":")

    def _get_interface_throats(self, phase1, phase2):
        conns = self.network.conns
        occ1 = self[f"pore.occupancy.{phase1}"][conns]
        occ2 = self[f"pore.occupancy.{phase2}"][conns]
        idx12, = np.where((occ1[:, 0] == 1) & (occ2[:, 1] == 1))
        idx21, = np.where((occ2[:, 0] == 1) & (occ1[:, 1] == 1))
        return idx12, idx21

    def _build_K(self):
        """Updates the global partition coefficient array"""
        prefix = self.settings["partition_coef_prefix"]
        self._K = np.ones(self.Nt, dtype=float)
        # Find all binary partition coefficient models
        models = [k for k in self.models.keys() if k.startswith(prefix)]
        # Modify the global partition coefficient for each phase pair
        for model in models:
            K12 = self[model]
            phase1, phase2 = self._get_phase_labels(model)
            idx12, idx21 = self._get_interface_throats(phase1, phase2)
            self._K[idx12] = K12[idx12]
            self._K[idx21] = 1 / K12[idx21]
        # Store a reference in self as a propname for convenience
        self[f"{prefix}.global"][:] = self._K

    def _interleave_data(self, prop):
        """Gathers property values from component phases to build a single array."""
        element = self._parse_element(prop)[0]
        vals = np.zeros(self._count(element=element), dtype=float)
        # Retrieve property from constituent phases (weight = occupancy)
        for phase in self.phases.values():
            vals += phase[prop] * self[f"{element}.occupancy.{phase.name}"]
        return vals

    def _set_automatic_throat_occupancy(self, mode="mean"):
        """
        Automatically interpolates throat occupancy based on that in
        adjacent pores. This method doesn't return any value.

        Parameters
        ----------
        mode : str
            Interpolation method. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'mean'       sets the throat occupancy as the average of that in
                         adjacent pores.
            'min'        sets the throat occupancy as the minimum value of
                         that in adjacent pores.
            'max'        sets the throat occupancy as the maximum value of
                         that in adjacent pores.
            ===========  =====================================================

        """
        self.settings['throat_occupancy'] = 'automatic'
        for phase in self.phases.values():
            self.add_model(propname=f"throat.occupancy.{phase.name}",
                           model=misc.from_neighbor_pores,
                           prop=f"pore.occupancy.{phase.name}",
                           mode=mode)


def multiphase_diffusion(phase,
                         pore_diffusivity="pore.diffusivity",
                         throat_diffusivity="throat.diffusivity",
                         size_factors="throat.diffusive_size_factors",
                         partition_coef_global="throat.partition_coef.global"):
    r"""
    Calculates the diffusive conductance of conduits for multiphase systems.

    Parameters
    ----------
    %(phase)s
    pore_diffusivity : str
        %(dict_blurb)s pore diffusivity
    throat_diffusivity : str
        %(dict_blurb)s throat diffusivity
    size_factors : str
        %(dict_blurb)s conduit size factors

    Returns
    -------
    %(return_arr)s diffusive conductance

    Notes
    -----
    This method assumes that ``phase["partition_coef"]`` contains information
    on binary phase partitioning. See ``MultiPhase`` class documentation for
    more information.

    """
    network = phase.network
    cn = network.conns
    SF = network[size_factors]
    if isinstance(SF, dict):
        F1, Ft, F2 = SF.values()
    elif SF.ndim > 1:
        F1, Ft, F2 = SF.T
    else:
        F1, Ft, F2 = np.inf, SF, np.inf

    # Fetch model parameters
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity]
    g1 = D1 * F1
    gt = Dt * Ft
    g2 = D2 * F2

    # Apply Henry's partitioning coefficient
    # Note: m12 = (G21*c1 - G12*c2)  NOT  (G12*c1 - G21*c2)
    K12 = phase[partition_coef_global]
    G21 = (1/g1 + 0.5/gt + K12 * (1/g2 + 0.5/gt)) ** -1
    G12 = K12 * G21

    return np.vstack((G12, G21)).T
