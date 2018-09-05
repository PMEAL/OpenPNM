import numpy as np
import scipy.sparse as sprs
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


class TransientReactiveTransport(ReactiveTransport):
    r"""
    A subclass of ReactiveTransport for transient/steady-state simulations

    Parameters
    ----------
    network : OpenPNM Network object
        The Network with which this algorithm is associated.

    project : OpenPNM Project object
        Either a Network or a Project must be specified.

    Notes
    -----

    This subclass performs steady and transient simulations of transport
    phenomena with reactions when source terms are added. It supports 3 time
    discretization schemes; 'steady' to perform a steady-state simulation, and
    'implicit' (fast, 1st order accurate) and 'cranknicolson' (slow, 2nd order
    accurate) both for transient simulations.
    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   't_initial': 0,
                   't_final': 10,
                   't_step': 0.1,
                   't_output': 1e+08,
                   't_tolerance': 1e-06,
                   'r_tolerance': 1e-04,
                   't_scheme': 'implicit',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_scheme': ''},
                           'set_IC':       {'values': None},
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
        self._A_steady = None  # Initialize the steady sys of eqs A matrix
        if phase is not None:
            self.setup(phase=phase)

    @docstr.get_sectionsf('TransientReactiveTransport.setup',
                          sections=['Parameters', 'Other Parameters'])
    @docstr.dedent
    def setup(self, phase=None, quantity='', conductance='',
              t_initial=None, t_final=None, t_step=None, t_output=None,
              t_tolerance=None, t_scheme='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings

        Parameters
        ----------
        %(ReactiveTransport.setup.parameters)s

        Other Parameters
        ----------------
        t_initial : scalar, smaller than 't_final'
            The simulation's start time. The default value is 0.
        t_final : scalar, bigger than 't_initial'
            The simulation's end time. The default value is 10.
        t_step : scalar, between 't_initial' and 't_final'
            The simulation's time step. The default value is 0.1.
        t_output : scalar
            Output interval to store transient solutions. The default value
            is 1e+08. Initial and steady-state (if reached) fields are always
            stored. If 't_output' > 't_final', no transient data is stored.
            If 't_output' is not a multiple of 't_step', 't_output' will be
            approximated.
        t_tolerance : scalar
            Transient solver tolerance. The simulation stops (before reaching
            't_final') when the residual falls below 't_tolerance'. The
            default value is 1e-06. The 'residual' measures the variation from
            one time-step to another in the value of the 'quantity' solved for.
        r_tolerance : scalar
            Tolerance to achieve within each time step. The solver passes to
            next time step when 'residual' falls below 'r_tolerance'. The
            default value is 1e-04.
        t_scheme : string
            The time discretization scheme. Three options available: 'steady'
            to perform a steady-state simulation, and 'implicit' (fast, 1st
            order accurate) and 'cranknicolson' (slow, 2nd order accurate) both
            for transient simulations. The default value is 'implicit'.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if t_initial:
            self.settings['t_initial'] = t_initial
        if t_final:
            self.settings['t_final'] = t_final
        if t_step:
            self.settings['t_step'] = t_step
        if t_output:
            self.settings['t_output'] = t_output
        if t_tolerance:
            self.settings['t_tolerance'] = t_tolerance
        if t_scheme:
            self.settings['t_scheme'] = t_scheme
        self.settings.update(kwargs)

    def set_IC(self, values):
        r"""
        A method to set simulation initial conditions

        Parameters
        ----------
        values : ND-array or scalar
            Set the initial conditions using an 'Np' long array. 'Np' being
            the number of pores. If a scalar is given, the same value is
            imposed to all pores.
        """
        self[self.settings['quantity']] = values
        converted_array = self[self.settings['quantity']].astype('float64')
        self[self.settings['quantity']] = converted_array

    def _t_update_A(self):
        r"""
        A method to update 'A' matrix at each time step according to 't_scheme'
        """
        network = self.project.network
        Vi = network['pore.volume']
        dt = self.settings['t_step']
        s = self.settings['t_scheme']
        if (s == 'implicit'):
            f1, f2 = 1, 1
        elif (s == 'cranknicolson'):
            f1, f2 = 0.5, 1
        elif (s == 'steady'):
            f1, f2 = 1, 0
        # Compute A (operations involve conversion to 'csr')
        A = ((f2/dt) * sprs.coo_matrix.multiply(
            sprs.coo_matrix(np.reshape(Vi, (self.Np, 1)), shape=(self.Np,)),
            sprs.identity(self.Np, format='coo')) + f1 * self._A_steady)
        # Convert A to 'coo' format to apply BCs
        A = sprs.coo_matrix(A)
        self._A = A
        return A

    def _t_update_b(self):
        r"""
        A method to update 'b' array at each time step according to
        't_scheme' and the source term value
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        Vi = network['pore.volume']
        dt = self.settings['t_step']
        s = self.settings['t_scheme']
        if (s == 'implicit'):
            f1, f2, f3 = 1, 1, 0
        elif (s == 'cranknicolson'):
            f1, f2, f3 = 0.5, 1, 0
        elif (s == 'steady'):
            f1, f2, f3 = 1, 0, 1
        x_old = self[self.settings['quantity']]
        b = (f2*(1-f1)*(-self._A_steady)*x_old +
             f2*(Vi/dt)*x_old +
             f3*np.zeros(shape=(self.Np, ), dtype=float))
        self._update_physics()
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Update b
            b[Ps] = b[Ps] - f2*(1-f1)*(phase[item+'.'+'rate'][Ps])
        self._b = b
        return b

    def run(self, t=None):
        r"""
        Builds 'A' matrix of the steady system of equations to be used at each
        time step to build transient 'A' and 'b'. Imposes the initial
        conditions and stores the initial field. Initialize transient 'A', 'b',
        and source term (if present) and finally calls the transient solver.

        Parameters
        ----------
        t : scalar
            The time to start the simulation from. If no time is specified, the
            simulation starts from 't_initial' defined in the settings.
        """
        logger.info('―'*80)
        logger.info('Running TransientTransport')
        # If solver used in steady mode, no need to add ICs
        if (self.settings['t_scheme'] == 'steady'):
            self[self.settings['quantity']] = 0.0
        # If ICs are not defined, show an error
        if self[self.settings['quantity']] is None:
            logger.error('Initial conditions not defined')
        # Save A matrix of the steady sys of eqs (WITHOUT BCs applied)
        self._A_steady = (self.A).copy()
        # Initialize A and b with BCs applied
        self._t_update_A()
        self._t_update_b()
        self._apply_BCs()
        self._A_t = (self._A).copy()
        self._b_t = (self._b).copy()
        if t is None:
            t = self.settings['t_initial']
        # Create S1 & S1 for 1st Picard's iteration
        self._update_physics()

        self._run_transient(t=t)

    def _run_transient(self, t):
        """r
        Performs a transient simulation according to the specified settings
        updating 'b' and calling '_t_run_reactive' at each time step.
        Stops after reaching the end time 't_final' or after achieving the
        specified tolerance 't_tolerance'. Stores the initial and steady-state
        (if obtained) fields in addition to transient data (according to the
        specified 't_output').

        Parameters
        ----------
        t : scalar
            The time to start the simulation from.

        Notes
        -----
        Transient solutions are stored on the object under
        ``pore.quantity_timeStepIndex`` where *quantity* is specified in the
        ``settings`` attribute. Initial field is stored as
        ``pore.quantity_initial``. Steady-state solution (if reached) is stored
        as ``pore.quantity_steady``. Current solution is stored as
        ``pore.quantity``.
        """
        tf = self.settings['t_final']
        dt = self.settings['t_step']
        to = self.settings['t_output']
        tol = self.settings['t_tolerance']
        s = self.settings['t_scheme']
        res_t = 1e+06  # Initialize the residual

        # Make sure 'tf' and 'to' are multiples of 'dt'
        tf = tf + (dt-(tf % dt))*((tf % dt) != 0)
        to = to + (dt-(to % dt))*((to % dt) != 0)
        self.settings['t_final'] = tf
        self.settings['t_output'] = to
        outputs = np.append(np.arange(t+to, tf, to), tf)

        if (s == 'steady'):  # If solver in steady mode, do one iteration
            logger.info('    Running in steady mode')
            x_old = self[self.settings['quantity']]
            self._t_run_reactive(x=x_old)
            x_new = self[self.settings['quantity']]

        else:  # Do time iterations
            # Export the initial field (t=t_initial)
            vals = self[self.settings['quantity']]
            self[self.settings['quantity']+'_initial'] = vals
            for time in np.arange(t+dt, tf+dt, dt):
                if (res_t >= tol):  # Check if the steady state is reached
                    logger.info('    Current time step: '+str(time)+' s')
                    x_old = self[self.settings['quantity']]

                    self._t_run_reactive(x=x_old)
                    x_new = self[self.settings['quantity']]

                    # Compute the residual
                    res_t = np.sum(np.absolute(x_old**2 - x_new**2))
                    logger.info('        Residual: '+str(res_t))
                    # Output transient solutions. Round time to ensure every
                    # value in outputs is exported.
                    if round(time, 12) in outputs:
                        ind = np.where(outputs == round(time, 12))[0][0]
                        self[self.settings['quantity']+'_'+str(ind)] = x_new
                        logger.info('        Exporting time step: ' +
                                    str(time)+' s')
                    # Update A and b and apply BCs
                    self._t_update_A()
                    self._t_update_b()
                    self._apply_BCs()
                    self._A_t = (self._A).copy()
                    self._b_t = (self._b).copy()

                else:  # Stop time iterations if residual < t_tolerance
                    self[self.settings['quantity'] + '_steady'] = x_new
                    logger.info('        Exporting time step: '+str(time)+' s')
                    break
            if (round(time, 12) == tf):
                logger.info('    Maximum time step reached: '+str(time)+' s')
            else:
                logger.info('    Transient solver converged after: ' +
                            str(time)+' s')

    def _t_run_reactive(self, x):
        """r
        Repeatedly updates transient 'A', 'b', and the solution guess within
        each time step according to the applied source term then calls '_solve'
        to solve the resulting system of linear equations. Stops when the
        residual falls below 'r_tolerance'.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        Returns
        -------
        x_new : ND-array
            Solution array.

        Notes
        -----
        Description of 'relaxation_quantity' and 'max_iter' settings can be
        found in the parent class 'ReactiveTransport' documentation.
        """
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        self[self.settings['quantity']] = x
        relax = self.settings['relaxation_quantity']
        res = 1e+06
        for itr in range(int(self.settings['max_iter'])):
            if res >= self.settings['r_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                self[self.settings['quantity']] = x
                self._A = (self._A_t).copy()
                self._b = (self._b_t).copy()
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
