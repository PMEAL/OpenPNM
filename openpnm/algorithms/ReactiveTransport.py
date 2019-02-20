import numpy as np
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
                   'r_tolerance': 0.001,
                   'max_iter': 5000,
                   'relaxation_source': 1,
                   'relaxation_quantity': 1,
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            'r_tolerance': None,
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

    def setup(self, phase=None, quantity='', conductance='', r_tolerance=None,
              max_iter=None, relaxation_source=None,
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

        r_tolerance : scalar
            Tolerance to achieve. The solver returns a solution when 'residual'
            falls below 'r_tolerance'. The default value is 0.001.

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
        if r_tolerance:
            self.settings['r_tolerance'] = r_tolerance
        if max_iter:
            self.settings['max_iter'] = max_iter
        if relaxation_source:
            self.settings['relaxation_source'] = relaxation_source
        if relaxation_quantity:
            self.settings['relaxation_quantity'] = relaxation_quantity
        super().setup(**kwargs)

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

        if (not np.all(np.isnan(self['pore.bc_value'][locs]))) or \
           (not np.all(np.isnan(self['pore.bc_rate'][locs]))):
            raise Exception('Boundary conditions already present in given ' +
                            'pores, cannot also assign source terms')
        self[propname] = locs
        self.settings['sources'].append(propname)

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
                raise Exception('Source term already present in given ' +
                                'pores, cannot also assign boundary ' +
                                'conditions')
        # Then call parent class function if above check passes
        super()._set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues,
                        mode=mode)

    _set_BC.__doc__ = GenericTransport._set_BC.__doc__

    def _update_physics(self):
        """r
        Update physics using the current value of 'quantity'

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the phase.
        This method was implemented relaxing one of the OpenPNM rules of
        algorithms not being able to write into phases.
        """
        phase = self.project.phases()[self.settings['phase']]
        physics = self.project.find_physics(phase=phase)
        for item in self.settings['sources']:
            # Regenerate models with new guess
            quantity = self.settings['quantity']
            # Put quantity on phase so physics finds it when regenerating
            phase[quantity] = self[quantity]
            # Regenerate models, on either phase or physics
            phase.regenerate_models(propnames=item)
            for phys in physics:
                phys.regenerate_models(propnames=item)

    def _apply_sources(self):
        """r
        Update 'A' and 'b' applying source terms to specified pores

        Notes
        -----
        Applying source terms to 'A' and 'b' is performed after (optionally)
        under-relaxing the source term to improve numerical stability. Physics
        are also updated before applying source terms to ensure that source
        terms values are associated with the current value of 'quantity'.
        In the case of a transient simulation, the updates in 'A' and 'b'
        also depend on the time scheme.
        """
        if self.settings['t_scheme'] == 'cranknicolson':
            f1 = 0.5
        else:
            f1 = 1
        phase = self.project.phases()[self.settings['phase']]
        relax = self.settings['relaxation_source']
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Add S1 to diagonal of A
            # TODO: We need this to NOT overwrite the A and b, but create
            # copy, otherwise we have to regenerate A and b on each loop
            datadiag = self._A.diagonal().copy()
            # Source term relaxation
            S1_old = phase[item+'.'+'S1'][Ps].copy()
            S2_old = phase[item+'.'+'S2'][Ps].copy()
            self._update_physics()
            S1 = phase[item+'.'+'S1'][Ps]
            S2 = phase[item+'.'+'S2'][Ps]
            S1 = relax*S1 + (1-relax)*S1_old
            S2 = relax*S2 + (1-relax)*S2_old
            phase[item+'.'+'S1'][Ps] = S1
            phase[item+'.'+'S2'][Ps] = S2
            datadiag[Ps] = datadiag[Ps] - f1*S1
            # Add S1 to A
            self._A.setdiag(datadiag)
            # Add S2 to b
            self._b[Ps] = self._b[Ps] + f1*S2

    def run(self, x=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        """
        logger.info('Running ReactiveTransport')
        # Create S1 & S1 for 1st Picard's iteration
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        self[self.settings['quantity']] = x
        self._update_physics()

        x = self._run_reactive(x=x)
        self[self.settings['quantity']] = x

    def _run_reactive(self, x):
        r"""
        Repeatedly updates 'A', 'b', and the solution guess within according
        to the applied source term then calls '_solve' to solve the resulting
        system of linear equations.
        Stops when the residual falls below 'r_tolerance' or when the maximum
        number of iterations is reached.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        Returns
        -------
        x_new : ND-array
            Solution array.
        """
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        self[self.settings['quantity']] = x
        relax = self.settings['relaxation_quantity']
        res = 1e+06  # Initialize the residual
        for itr in range(int(self.settings['max_iter'])):
            if res >= self.settings['r_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                self[self.settings['quantity']] = x
                self._build_A(force=True)
                self._build_b(force=True)
                self._apply_BCs()
                self._apply_sources()
                x_new = self._solve()
                # Relaxation
                x_new = relax*x_new + (1-relax)*self[self.settings['quantity']]
                self[self.settings['quantity']] = x_new
                res = np.sum(np.absolute(x**2 - x_new**2))
                x = x_new
            if (res < self.settings['r_tolerance'] or
                    self.settings['sources'] == []):
                logger.info('Solution converged: ' + str(res))
                break
        return x_new
