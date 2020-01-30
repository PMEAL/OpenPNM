import numpy as np
from numpy.linalg import norm
from openpnm.algorithms import GenericTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class ReactiveTransport(GenericTransport):
    r"""
    A subclass for steady-state simulations with (optionally) source terms

    Parameters
    ----------
    network : OpenPNM Network object
        The Network with which this algorithm is associated.

    project : OpenPNM Project object
        Either a Network or a Project must be specified.

    Notes
    -----

    This subclass performs steady simulations of transport phenomena with
    reactions when source terms are added.
    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'sources': [],
                   'max_iter': 5000,
                   'relaxation_source': 1.0,
                   'relaxation_quantity': 1.0,
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            'rxn_tolerance': None,
                                            'max_iter': None,
                                            'relaxation_source': None,
                                            'relaxation_quantity': None},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='',
              rxn_tolerance=None, max_iter=None, relaxation_source=None,
              relaxation_quantity=None, **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run. If no value is
            given, the existing value is kept.

        quantity : string
            The name of the physical quantity to be calcualted such as
            ``'pore.xxx'``.

        conductance : string
            The name of the pore-scale transport conductance values. These
            are typically calculated by a model attached to a *Physics* object
            associated with the given *Phase*. Example; ``'throat.yyy'``.

        rxn_tolerance : scalar
            Tolerance to achieve. The solver returns a solution when 'residual'
            falls below 'rxn_tolerance'. The default value is 1e-05.

        max_iter : scalar
            The maximum number of iterations the solver can perform to find
            a solution. The default value is 5000.

        relaxation_source : scalar, between 0 and 1
            A relaxation factor to control under-relaxation of the source term.
            Factor approaching 0 : improved stability but slow simulation.
            Factor approaching 1 : fast simulation but may be unstable.
            Default value is 1 (no under-relaxation).

        relaxation_quantity :  scalar, between 0 and 1
            A relaxation factor to control under-relaxation for the quantity
            solving for.
            Factor approaching 0 : improved stability but slow simulation.
            Factor approaching 1 : fast simulation but may be unstable.
            Default value is 1 (no under-relaxation).

        Notes
        -----
        Under-relaxation is a technique used for improving stability of a
        computation, particularly in the presence of highly non-linear terms.
        Under-relaxation used here limits the change in a variable from one
        iteration to the next. An optimum choice of the relaxation factor is
        one that is small enough to ensure stable simulation and large enough
        to speed up the computation.
        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if rxn_tolerance:
            self.settings['rxn_tolerance'] = rxn_tolerance
        if max_iter:
            self.settings['max_iter'] = max_iter
        if relaxation_source:
            self.settings['relaxation_source'] = relaxation_source
        if relaxation_quantity:
            self.settings['relaxation_quantity'] = relaxation_quantity
        super().setup(**kwargs)

    def reset(self, source_terms=True, **kwargs):
        r"""
        """
        super().reset(**kwargs)
        if source_terms:
            for item in self.settings['sources']:
                self.pop(item)
            self.settings.pop('sources', None)

    def set_source(self, propname, pores):
        r"""
        Applies a given source term to the specified pores

        Parameters
        ----------
        propname : string
            The property name of the source term model to be applied

        pores : array_like
            The pore indices where the source term should be applied

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
        # Set source term
        self[propname] = locs
        self.settings['sources'].append(propname)
        # Add source term as an iterative prop
        self.set_iterative_props(propname)

    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        Apply boundary conditions to specified pores if no source terms are
        already assigned to these pores. Otherwise, raise an error.

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied

        bctype : string
            Specifies the type or the name of boundary condition to apply. The
            types can be one one of the following:

            - *'value'* : Specify the value of the quantity in each location
            - *'rate'* : Specify the flow rate into each location

        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.  If
            a single value is given, it's assumed to apply to all locations.
            Different values can be applied to all pores in the form of an
            array of the same length as ``pores``.

        mode : string, optional
            Controls how the conditions are applied.  Options are:

            *'merge'*: (Default) Adds supplied boundary conditions to already
            existing conditions.

            *'overwrite'*: Deletes all boundary condition on object then add
            the given ones

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm. Use ``remove_BCs`` to
        clear existing BCs before applying new ones or ``mode='overwrite'``
        which removes all existing BC's before applying the new ones.
        """
        # First check that given pores do not have source terms already set
        for item in self.settings['sources']:
            if np.any(self[item][pores]):
                raise Exception('Source term already present in given '
                                + 'pores, cannot also assign boundary '
                                + 'conditions')
        # Then call parent class function if above check passes
        super()._set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues,
                        mode=mode)

    _set_BC.__doc__ = GenericTransport._set_BC.__doc__

    def _update_iterative_props(self):
        """r
        Update physics using the current value of ``quantity``

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the phase.
        This method was implemented relaxing one of the OpenPNM rules of
        algorithms not being able to write into phases.
        """
        phase = self.project.phases()[self.settings['phase']]
        physics = self.project.find_physics(phase=phase)
        # Put quantity on phase so physics finds it when regenerating
        phase.update(self.results())
        # Regenerate iterative props with new guess
        iterative_props = self.settings["iterative_props"]
        phase.regenerate_models(propnames=iterative_props)
        for physic in physics:
            physic.regenerate_models(iterative_props)

    def _apply_sources(self):
        """r
        Update ``A`` and ``b`` applying source terms to specified pores

        Notes
        -----
        Applying source terms to ``A`` and ``b`` is performed after (optionally)
        under-relaxing the source term to improve numerical stability. Physics
        are also updated before applying source terms to ensure that source
        terms values are associated with the current value of 'quantity'.

        Warnings
        --------
        In the case of a transient simulation, the updates in ``A`` and ``b``
        also depend on the time scheme. So, ``_correct_apply_sources()`` needs to
        be run afterwards to correct the already applied relaxed source terms.
        """
        phase = self.project.phases()[self.settings['phase']]
        w = self.settings['relaxation_source']

        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Source term relaxation
            S1, S2 = [phase[item + '.' + x][Ps] for x in ['S1', 'S2']]
            # Get old values of S1 and S2
            try:
                X1, X2 = [phase[item + '.' + x + '.old'][Ps] for x in ['S1', 'S2']]
            # S1.old and S2.old are not yet available in 1st iteration
            except KeyError:
                X1, X2 = S1.copy(), S2.copy()
            S1 = phase[item + '.' + 'S1'][Ps] = w * S1 + (1-w) * X1
            S2 = phase[item + '.' + 'S2'][Ps] = w * S2 + (1-w) * X2
            # Add "relaxed" S1 and S2 to A and b
            datadiag = self._A.diagonal().copy()
            datadiag[Ps] = datadiag[Ps] - S1
            self._A.setdiag(datadiag)
            self._b[Ps] = self._b[Ps] + S2

        # Replace old values of S1 and S2 by their current values
        for item in self.settings['sources']:
            phase[item + '.' + 'S1.old'] = phase[item + '.' + 'S1'].copy()
            phase[item + '.' + 'S2.old'] = phase[item + '.' + 'S2'].copy()

    def run(self, x0=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable
        """
        quantity = self.settings['quantity']
        logger.info('Running ReactiveTransport')

        # Create S1 & S2 for the 1st Picard iteration
        if x0 is None:
            x0 = np.zeros(self.Np, dtype=float)
        # Write initial guess to algorithm obj (for _update_iterative_props to work)
        self[quantity] = x0
        x = self._run_reactive(x0)
        self[quantity] = x

    def _run_reactive(self, x0):
        r"""
        Repeatedly updates ``A``, ``b``, and the solution guess within according
        to the applied source term then calls ``_solve`` to solve the resulting
        system of linear equations.

        Stops when the residual falls below ``solver_tol * norm(b)`` or when
        the maximum number of iterations is reached.

        Parameters
        ----------
        x0 : ND-array
            Initial guess of unknown variable

        Returns
        -------
        x : ND-array
            Solution array.
        """
        x = x0
        w = self.settings['relaxation_quantity']
        quantity = self.settings['quantity']
        max_it = self.settings['max_iter']

        for itr in range(max_it):
            # Update iterative properties on phase and physics
            self._update_iterative_props()
            # Build A and b, apply BCs/source terms
            self._build_A()
            self._build_b()
            self._apply_BCs()
            self._apply_sources()
            # Check solution convergence
            res = self._get_residual()
            if self._is_converged():
                logger.info(f'Solution converged: {res:.4e}')
                return x
            logger.info(f'Tolerance not met: {res:.4e}')
            # Solve, use relaxation, and update solution on algorithm obj
            self[quantity] = x = self._solve(x0=x) * w + x * (1 - w)

        # Check solution convergence after max_it iterations
        if not self._is_converged():
            raise Exception("Not converged after {max_it} iterations.")

    def _is_converged(self):
        r"""
        Check if solution has converged based on the following criterion:
            res <= max(norm(b) * tol, atol)
        """
        res = self._get_residual()
        # Verify that residual is finite (i.e. not inf/nan)
        if not np.isfinite(res):
            logger.warning(f'Solution diverged: {res:.4e}')
            raise Exception(f"Solution diverged, undefined residual: {res:.4e}")
        # Check convergence
        tol = self.settings["solver_tol"]
        res_tol = norm(self.b) * tol
        flag_converged = True if res <= res_tol else False
        return flag_converged
