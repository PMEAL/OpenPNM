import numpy as np
from numpy.linalg import norm
from scipy.optimize.nonlin import TerminationCondition
# Uncomment this line when we stop supporting Python 3.6
# from dataclasses import dataclass, field
# from typing import List
from openpnm.algorithms import GenericTransport
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


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
    relaxation_quantity : float (default = 1.0)
        A relaxation factor to control under-relaxation for the quantity
        solving for. Factor approaching 0 leads to improved stability but
        slower simulation. Factor approaching 1 gives fast simulation but
        may be unstable.
    newton_maxiter : int
        Maximum number of iterations allowed for the nonlinear solver to
        converge.

    ----

    **The following parameters pertain to the ``GenericTransport`` class**

    %(GenericTransportSettings.other_parameters)s

    """

    newton_maxiter = 5000
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

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(ReactiveTransportSettings)
        self.settings.update(settings)

    @docstr.dedent
    def reset(self, source_terms=False, **kwargs):
        r"""
        %(GenericTransport.reset.full_desc)s

        Parameters
        ----------
        %(GenericTransport.reset.parameters)s
        source_terms : bool
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

    def set_source(self, propname, pores, mode='overwrite'):
        r"""
        Applies a given source term to the specified pores.

        Parameters
        ----------
        propname : str
            The property name of the source term model to be applied.
        pores : array_like
            The pore indices where the source term should be applied.
        mode : str
            Controls how the sources are applied (see table under Notes).
            The default is 'overwrite'. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'merge'      Adds supplied source term to already existing
                         ones.
            'overwrite'  Deletes all existing source terms of the given
                         ``propname`` then adds the specified new ones.
            ===========  =====================================================

        Notes
        -----
        Source terms cannot be applied in pores where boundary conditions
        have already been set. Attempting to do so will result in an error
        being raised.

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
        Phase and physics objects are also updated before applying source
        terms to ensure that source terms values are associated with the
        current value of 'quantity'.

        """
        phase = self.project.find_phase(self)
        for item in self.settings['sources']:
            # Fetch linearized values of the source term
            Ps = self.pores(item)
            S1, S2 = [phase[f"{item}.{Si}"] for Si in ["S1", "S2"]]
            # Modify A and b: diag(A) += -S1, b += S2
            diag = self.A.diagonal()
            diag[Ps] += -S1[Ps]
            self.A.setdiag(diag)
            self.b[Ps] += S2[Ps]

    def _run_special(self, solver, x0):
        r"""
        Repeatedly updates ``A``, ``b``, and the solution guess within
        according to the applied source term then calls ``_solve`` to
        solve the resulting system of linear equations.

        Stops when the max-norm of the residual drops by at least
        ``f_rtol``:

            ``norm(R_n) < norm(R_0) * f_rtol``

        AND

            ``norm(dx) < norm(x) * x_rtol``

        where R_i is the residual at ith iteration, x is the solution at
        current iteration, and dx is the change in the solution between two
        consecutive iterations. ``f_rtol`` and ``x_rtol`` are defined in
        the algorithm's settings under: ``alg.settings['f_rtol']``, and
        ``alg.settings['x_rtol']``, respectively.

        Parameters
        ----------
        x0 : ndarray
            Initial guess of the unknown variable

        """
        w = self.settings['relaxation_quantity']
        quantity = self.settings['quantity']
        maxiter = self.settings['newton_maxiter']
        f_rtol = self.settings['f_rtol']
        x_rtol = self.settings['f_xtol']
        xold = self[quantity]
        condition = TerminationCondition(f_rtol=f_rtol, x_rtol=x_rtol)

        for i in range(maxiter):
            dx = self[quantity] - xold
            res = self._get_residual()
            is_converged = condition.check(f=res, x=xold, dx=dx)
            if is_converged:
                logger.info(f'Solution converged, residual norm: {norm(res):.4e}')
                return
            super()._run_special(solver=solver, x0=xold, w=w)
            xold = self[quantity]
            logger.info(f'Iteration #{i:<4d} | Residual norm: {norm(res):.4e}')
        logger.critical(f"{self.name} didn't converge after {maxiter} iterations")

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

    def _get_residual(self, x=None):
        r"""
        Calculate solution residual based on the given ``x`` based on the
        following formula:
            ``res = norm(A*x - b)``
        """
        if x is None:
            quantity = self.settings['quantity']
            x = self[quantity]
        return self.A * x - self.b

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
