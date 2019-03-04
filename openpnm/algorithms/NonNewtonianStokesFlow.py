import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class NonNewtonianStokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the
    network.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.pressure',
                   'conductance': 'throat.nonNewtonian_hydraulic_conductance',
                   'tolerance': 1e-5,
                   'max_iter': 10,
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': ''},
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

    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run.  If no value is
            given, the existing value is kept.

        quantity : string
            The name of the physical quantity to be calcualted.  If no value is
            given, the existing value is kept.  The default value is
            ``'pore.pressure'``.

        conductance : string
            The name of the pore-scale transport conductance values.  These
            are typically calculate by a model attached to a *Physics* object
            associated with the given *Phase*.  If no value is given, the
            existing value is kept.  The default value is
            ``'throat.hydraulic_conductance'``.

        Notes
        -----
        Any additional arguments are added to the ``settings`` dictionary of
        the object.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        super().setup(**kwargs)

    def run(self):
        phase = self.project.phases()[self.settings['phase']]
        phys = self.project.find_physics(phase=phase)

        # Define initial conditions (if not defined by the user)
        try:
            self[self.settings['quantity']]
        except KeyError:
            self[self.settings['quantity']] = np.zeros(shape=[self.Np, ],
                                                       dtype=float)

        # Define tolerance and initialize residuals
        tol = self.settings['tolerance']
        res = 1e+06

        # Iterate until solution converges
        for itr in range(int(self.settings['max_iter'])):
            logger.info('Iter: ' + str(itr) + ', Res: ' + str(res))
            convergence = res < tol
            if not convergence:
                self._update_physics()
                phys[0].regenerate_models()
                p_old = self[self.settings['quantity']].copy()
                self._run_reactive(x=p_old)
                p_new = self[self.settings['quantity']].copy()
                # Residual
                res = np.sum(np.absolute(p_old**2 - p_new**2))
                phase.update(self.results())

            if convergence:
                logger.info('Solution converged: ' + str(res))
                break
