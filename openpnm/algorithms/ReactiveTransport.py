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
        if not source_terms:
            return
        # Remove item from label dictionary
        for item in self.settings['sources']:
            self.pop(item)
        # Reset the settings dict
        self.settings['sources'] = []
        # TODO: clean up S1/S2.old temp variables stored on self

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
            Controls how the sources are applied (see table under Notes).
            The default is 'overwrite'.

        Notes
        -----
        The following ``mode`` values are supported:

        ===========  =====================================================
        mode         meaning
        ===========  =====================================================
        'merge'      Adds supplied source term to already existing ones.
        'overwrite'  Deletes all existing source terms of the given
                     ``propname`` then adds the specified new ones.
        ===========  =====================================================

        Source terms cannot be applied in pores where boundary conditions have
        already been set. Attempting to do so will result in an error being
        raised.

        """
        propname = self._parse_prop(propname, "pore")
        locs = self.tomask(pores=pores)
        # Check if any BC is already set in the same locations
        locs_BC = np.isfinite(self['pore.bc_value']) + np.isfinite(self['pore.bc_rate'])
        if (locs & locs_BC).any():
            raise Exception("BCs present in given pores, can't assign source term")
        self.set_label(propname, pores=locs, mode=mode)
        # Check if propname already in source term list
        if propname not in self.settings['sources']:
            self.settings['sources'].append(propname)
        # Initialize '_propname.S1/S2.old' for optional under-relaxation
        _propname = "._".join(propname.split("."))
        self[f"{_propname}.S1.old"] = self[f"{_propname}.S2.old"] = 0.0

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
        locs = self.tomask(pores=pores or self.Ps)
        self.set_label(propname, pores=locs, mode='remove')
        # TODO: if pores=None: remove the label -> reuse in reset method

    def _update_iterative_props(self):
        """r
        Regenerates phase, geometries, and physics objects using the
        current value of ``quantity``.

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the
        phase, which is against one of the OpenPNM rules of objects not
        being able to write into each other.

        """
        iterative_props = self._get_iterative_props()
        if not iterative_props:
            return
        # Fetch objects associated with the algorithm
        phase = self.project.find_phase(self)
        physics = self.project.find_physics(phase=phase)
        geometries = self.project.geometries().values()
        # Update 'quantity' on phase with the most recent value
        quantity = self.settings['quantity']
        phase[quantity] = self[quantity]
        # Regenerate all associated objects
        phase.regenerate_models(propnames=iterative_props)
        for geom in geometries:
            geom.regenerate_models(iterative_props)
        for phys in physics:
            phys.regenerate_models(iterative_props)

    def _apply_sources(self):
        """r
        Updates ``A`` and ``b``, applying source terms to specified pores.

        Notes
        -----
        Applying source terms to ``A`` and ``b`` is performed after
        (optionally) under-relaxing the source term to improve numerical
        stability.

        Physics are also updated before applying source terms
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
        In the case of a transient simulation, the updates in ``A``, ``b``
        also depend on the time scheme. So, ``_correct_apply_sources()``
        needs to be run afterwards to correct the already applied relaxed
        source terms.

        """
        phase = self.project.find_phase(self)
        w = self.settings['relaxation_source']
        for item in self.settings['sources']:
            element, prop = item.split(".")
            _item = ".".join([element, "_" + prop])
            Ps = self.pores(item)
            S1, S2 = [phase[f"{item}.{Si}"] for Si in ["S1", "S2"]]
            X1, X2 = [self[f"{_item}.{Xi}"] for Xi in ["S1.old", "S2.old"]]
            # Source term relaxation
            S1 = phase[f"{item}.S1"][:] = w * S1 + (1 - w) * X1
            S2 = phase[f"{item}.S2"][:] = w * S2 + (1 - w) * X2
            # Modify A and b based on "relaxed" S1/S2
            diag = self.A.diagonal()
            diag[Ps] += -S1[Ps]
            self.A.setdiag(diag)
            self.b[Ps] += S2[Ps]
            # Replace old values of S1/S2 by their current values
            X1[:], X2[:] = S1, S2

    def _run_special(self, solver, x0):
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

        """
        w = self.settings['relaxation_quantity']
        quantity = self.settings['quantity']
        max_iter = self.settings['nlin_max_iter']
        for itr in range(max_iter):
            super()._run_special(solver=solver, x0=x0, w=w)
            x0 = self[quantity]  # Update x0 for next iteration
            if self._is_converged():
                logger.info(f'Solution converged: {self._get_residual():.4e}')
                return
            logger.info(f'Tolerance not met: {self._get_residual():.4e}')
        raise
        logger.critical(f"Not converged after {max_iter} iterations.")

    def _update_A_and_b(self):
        r"""
        Builds/updates A, b based on the recent solution on algorithm object.
        """
        self._update_iterative_props()
        super()._update_A_and_b()
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
        phase = self.project.find_phase(self)
        physics = self.project.find_physics(phase=phase)
        geometries = self.project.geometries().values()
        # Generate global dependency graph
        dg = nx.compose_all([x.models.dependency_graph(deep=True)
                             for x in [phase, *geometries, *physics]])
        base_props = [self.settings["quantity"]]
        # Find all props downstream that depend on 'quantity'
        dg = nx.DiGraph(nx.edge_dfs(dg, source=base_props))
        if len(dg.nodes) == 0:
            return []
        iterative_props = list(nx.dag.lexicographical_topological_sort(dg))
        # Remove 'quantity' from iterative_props since it is not!
        iterative_props.remove(self.settings["quantity"])
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
        msg = "Source term already present in given pores, can't assign BCs"
        # Ensure that given pores do not have source terms already set
        for item in self.settings['sources']:
            if np.any(self[item][pores]):
                raise Exception(msg)
        # Assign BCs if above check passes
        super()._set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues, mode=mode)
