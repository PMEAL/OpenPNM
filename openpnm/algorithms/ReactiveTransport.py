import numpy as np
from openpnm.algorithms import GenericTransport
# Uncomment this line when we stop supporting Python 3.6
# from dataclasses import dataclass, field
# from typing import List
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


# class RelaxationSettings(GenericSettings):
#     r"""
#     This class is a demonstration of how we can add nested settings classes
#     to other settings classes to make categories for some settings.  This is
#     being appended to the ReactiveTransportSettings class under the
#     'relaxation' attribute, and it works as planned by allowing the nested
#     dot access to its parameters. More work would be required to get it
#     functional such as dealing with deeply nested dicts and so on, but it
#     works in principal.
#     """
#     source = 1.0
#     quantity = 1.0


@docstr.get_sections(base='ReactiveTransportSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
# Uncomment this line when we stop supporting Python 3.6
# @dataclass
class ReactiveTransportSettings(GenericSettings):
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s

    quantity : str
        The name of the physical quantity to be calculated
    conductance : str
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------
    sources : list
        List of source terms that have been added
    relaxation_source : float (default = 1.0)
        A relaxation factor to control under-relaxation of the source term.
        Factor approaching 0 leads to improved stability but slower simulation.
        Factor approaching 1 gives fast simulation but may be unstable.
    relaxation_quantity : float (default = 1.0)
        A relaxation factor to control under-relaxation for the quantity
        solving for. Factor approaching 0 leads to improved stability but
        slower simulation. Factor approaching 1 gives fast simulation but
        may be unstable.
    nlin_max_iter : int
        Maximum number of iterations allowed for the nonlinear solver to
        converge. This parameter is different that ``GenericTransport``'s
        ``solver_max_iter``.

    ----

    **The following parameters pertain to the ``GenericTransport`` class**

    %(GenericTransportSettings.other_parameters)s

    """

    nlin_max_iter = 5000
    # relaxation = RelaxationSettings()
    relaxation_source = 1.0
    relaxation_quantity = 1.0
    # Swap the following 2 lines when we stop supporting Python 3.6
    # sources: List = field(default_factory=lambda: [])
    sources = []


@docstr.get_sections(base='ReactiveTransport', sections=['Parameters'])
@docstr.dedent
class ReactiveTransport(GenericTransport):
    r"""
    A subclass for steady-state simulations with (optional) source terms.

    Parameters
    ----------
    %(GenericTransport.parameters)s

    Notes
    -----
    This subclass performs steady simulations of transport phenomena with
    reactions when source terms are added.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(ReactiveTransportSettings)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    @docstr.get_sections(base='ReactiveTransport.setup',
                         sections=['Parameters', 'Notes'])
    @docstr.dedent
    def setup(self, phase=None, quantity='', conductance='',
              nlin_max_iter=None, relaxation_source=None,
              relaxation_quantity=None, **kwargs):
        r"""
        This method takes several arguments that are essential to running
        the algorithm and adds them to the settings.

        Parameters
        ----------
        %(GenericTransportSettings.parameters)s
        %(ReactiveTransportSettings.parameters)s

        Notes
        -----
        Under-relaxation is a technique used for improving stability of a
        computation, particularly in the presence of highly non-linear
        terms. Under-relaxation used here limits the change in a variable
        from one iteration to the next. An optimum choice of the
        relaxation factor is one that is small enough to ensure stable
        simulation and large enough to speed up the computation.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if nlin_max_iter:
            self.settings['nlin_max_iter'] = nlin_max_iter
        if relaxation_source:
            self.settings['relaxation_source'] = relaxation_source
        if relaxation_quantity:
            self.settings['relaxation_quantity'] = relaxation_quantity
        super().setup(**kwargs)

    def run(self, x0=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x0 : ND-array
            Initial guess of unknown variable
        """
        self._validate_settings()
        # Check if A and b are well-defined
        self._validate_data_health()
        quantity = self.settings['quantity']
        logger.info('Running ReactiveTransport')
        x0 = np.zeros(self.Np, dtype=float) if x0 is None else x0
        self["pore.initial_guess"] = x0
        x = self._run_reactive(x0)
        self[quantity] = x

    @docstr.dedent
    def reset(self, source_terms=False, **kwargs):
        r"""
        %(GenericTransport.reset.full_desc)s

        Parameters
        ----------
        %(GenericTransport.reset.parameters)s
        source_terms : boolean
            If ``True`` removes source terms. The default is ``False``.

        """
        super().reset(**kwargs)
        if source_terms:
            # Remove item from label dictionary
            for item in self.settings['sources']:
                self.pop(item)
            # Reset the settings dict
            self.settings['sources'] = []

    def set_source(self, propname, pores, mode='overwrite'):
        r"""
        Applies a given source term to the specified pores.

        Parameters
        ----------
        propname : string
            The property name of the source term model to be applied.
        pores : array_like
            The pore indices where the source term should be applied.
        mode : str
            Controls how the sources are applied. Options are:

            'merge' - Adds supplied source term to already existing ones.
            'overwrite' - (default) Deletes all existing source terms of the
            given ``propname`` then adds the specified new ones

        Notes
        -----
        Source terms cannot be applied in pores where boundary conditions have
        already been set. Attempting to do so will result in an error being
        raised.

        """
        locs = self.tomask(pores=pores)
        # Check if any BC is already set in the same locations
        locs_BC = np.isfinite(self['pore.bc_value']) + np.isfinite(self['pore.bc_rate'])
        if (locs & locs_BC).any():
            raise Exception('Boundary conditions already present in given '
                            + 'pores, cannot also assign source terms')
        if propname not in self.keys():
            self[propname] = False
        if mode == 'merge':
            self[propname][locs] = True
        elif mode == 'overwrite':
            self[propname] = locs
        # Check if propname already in source term list
        if propname not in self.settings['sources']:
            self.settings['sources'].append(propname)

    def remove_source(self, propname, pores=None):
        r"""
        Removes source terms from specified pores.

        Parameters
        ----------
        propname : str
            The property name of the source term model to be removed.
        pores : array_like
            The pore indices where the source term should be applied.

        """
        if pores is None:
            pores = self.Ps
        locs = self.tomask(pores=pores)
        if propname not in self.keys():
            self[propname] = False
        self[propname][locs] = False

    def _update_iterative_props(self):
        """r
        Update physics using the current value of ``quantity``.

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the
        phase. This method was implemented relaxing one of the OpenPNM
        rules of algorithms not being able to write into phases.

        """
        phase = self.project.phases()[self.settings['phase']]
        physics = self.project.find_physics(phase=phase)
        geometries = self.project.geometries().values()
        # Regenerate iterative props with new guess
        iterative_props = self._get_iterative_props()
        if len(iterative_props) > 0:
            # Put quantity on phase so physics finds it when regenerating
            key = self.settings['quantity']
            phase[key] = self[key]
            phase.regenerate_models(propnames=iterative_props)
            for geometry in geometries:
                geometry.regenerate_models(iterative_props)
            for phys in physics:
                phys.regenerate_models(iterative_props)

    def _apply_sources(self):
        """r
        Update ``A`` and ``b`` applying source terms to specified pores.

        Notes
        -----
        Applying source terms to ``A`` and ``b`` is performed after
        (optionally) under-relaxing the source term to improve numerical
        stability. Physics are also updated before applying source terms
        to ensure that source terms values are associated with the current
        value of 'quantity'.

        For source term under-relaxation, old values of S1 and S2 need
        to be stored somewhere, we chose to store them on the algorithm
        object. This is because storing them on phase/physics creates
        unintended problems, ex. storing them on physics -> IO complains
        added depth to the NestedDict, and storing them on the phase
        object results in NaNs in case source term is only added to a
        subset of nodes, which breaks our _check_for_nans algorithm.

        Warnings
        --------
        In the case of a transient simulation, the updates in ``A`` and ``b``
        also depend on the time scheme. So, ``_correct_apply_sources()`` needs to
        be run afterwards to correct the already applied relaxed source terms.

        """
        phase = self.project.phases()[self.settings['phase']]
        w = self.settings['relaxation_source']

        for item in self.settings['sources']:
            element, prop = item.split(".")
            _item = ".".join([element, "_" + prop])
            first_iter = False if _item + ".S1.old" in self.keys() else True
            Ps = self.pores(item)
            # Fetch S1/S2 and their old values (don't exist on 1st iter)
            S1 = phase[item + ".S1"][Ps]
            S2 = phase[item + ".S2"][Ps]
            X1 = self[_item + ".S1.old"][Ps] if not first_iter else S1
            X2 = self[_item + ".S2.old"][Ps] if not first_iter else S2
            # Source term relaxation
            S1 = phase[item + '.S1'][Ps] = w * S1 + (1.0 - w) * X1
            S2 = phase[item + '.S2'][Ps] = w * S2 + (1.0 - w) * X2
            # Modify A and b based on "relaxed" S1/S2
            datadiag = self._A.diagonal().copy()
            datadiag[Ps] = datadiag[Ps] - S1
            self._A.setdiag(datadiag)
            self._b[Ps] = self._b[Ps] + S2
            # Replace old values of S1/S2 by their current values
            self[_item + ".S1.old"] = phase[item + ".S1"]
            self[_item + ".S2.old"] = phase[item + ".S2"]

    def _run_reactive(self, x0):
        r"""
        Repeatedly updates ``A``, ``b``, and the solution guess within
        according to the applied source term then calls ``_solve`` to
        solve the resulting system of linear equations.

        Stops when the residual falls below ``solver_tol * norm(b)`` or
        when the maximum number of iterations is reached.

        Parameters
        ----------
        x0 : ndarray
            Initial guess of unknown variable

        Returns
        -------
        x : ndarray
            Solution array.

        Notes
        -----
        The algorithm must at least complete one iteration, and hence the
        check for itr >= 1, because otherwise, _check_for_nans() never
        gets called in case there's something wrong with the data, and
        therefore, the user won't get notified about the root cause of the
        algorithm divergence.

        """
        w = self.settings['relaxation_quantity']
        quantity = self.settings['quantity']
        max_it = self.settings['nlin_max_iter']
        # Write initial guess to algorithm obj (for _update_iterative_props to work)
        self[quantity] = x = x0
        # Update A and b based on self[quantity]
        self._update_A_and_b()
        # Just in case you got a lucky guess, i.e. x0!
        if self._is_converged():
            logger.info(f'Solution converged: {self._get_residual():.4e}')
            return x

        for itr in range(max_it):
            # Solve, use relaxation, and update solution on algorithm obj
            self[quantity] = x = self._solve(x0=x) * w + x * (1 - w)
            self._update_A_and_b()
            # Check solution convergence
            if self._is_converged():
                logger.info(f'Solution converged: {self._get_residual():.4e}')
                return x
            logger.info(f'Tolerance not met: {self._get_residual():.4e}')

        if not self._is_converged():
            raise Exception(f"Not converged after {max_it} iterations.")

    def _update_A_and_b(self):
        r"""
        Updates A and b based on the most recent solution stored on
        algorithm object.
        """
        # Update iterative properties on phase, geometries, and physics
        self._update_iterative_props()
        # Build A and b, apply BCs/source terms
        self._build_A()
        self._build_b()
        self._apply_BCs()
        self._apply_sources()

    def _get_iterative_props(self):
        r"""
        Find and return properties that need to be iterated while running
        the algorithm.

        Notes
        -----
        This method was moved from ReactiveTransport class to
        GenericTransport because source terms are not necessarily the only
        properties that need iteration during an algorithm (ex.
        concentration-dependent conductance)

        """
        import networkx as nx
        phase = self.project.phases(self.settings['phase'])
        physics = self.project.find_physics(phase=phase)
        geometries = self.project.geometries().values()
        # Combine dependency graphs of phase and all physics/geometries
        dg = phase.models.dependency_graph(deep=True)
        for g in geometries:
            dg = nx.compose(dg, g.models.dependency_graph(deep=True))
        for p in physics:
            dg = nx.compose(dg, p.models.dependency_graph(deep=True))
        base_props = [self.settings["quantity"]]
        if base_props is None:
            return []
        # Find all props downstream that rely on "quantity"
        dg = nx.DiGraph(nx.edge_dfs(dg, source=base_props))
        if len(dg.nodes) == 0:
            return []
        iterative_props = list(nx.dag.lexicographical_topological_sort(dg))
        # "quantity" shouldn't be in the returned list but "variable_props" should
        try:
            iterative_props.remove(self.settings["quantity"])
        except ValueError:
            pass
        return iterative_props

    @docstr.dedent
    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        Apply boundary conditions to specified pores if no source terms
        are already assigned to these pores. Otherwise, raise an error.

        Parameters
        ----------
        %(GenericTransport._set_BC.parameters)s

        Notes
        -----
        %(GenericTransport._set_BC.notes)s

        """
        # First check that given pores do not have source terms already set
        for item in self.settings['sources']:
            if np.any(self[item][pores]):
                raise Exception(
                    'Source term already present in given pores, cannot also'
                    ' assign boundary conditions'
                )
        # Then call parent class function if above check passes
        super()._set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues, mode=mode)
