import numpy as np
import scipy.sparse as sprs
from decimal import Decimal as dc
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


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
                   'rxn_tolerance': 1e-05,
                   't_precision': 12,
                   't_scheme': 'implicit',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_precision': None,
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

    def setup(self, phase=None, quantity='', conductance='',
              t_initial=None, t_final=None, t_step=None, t_output=None,
              t_tolerance=None, t_precision=None, t_scheme='', **kwargs):
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

        t_initial : scalar, smaller than 't_final'
            The simulation's start time. The default value is 0.

        t_final : scalar, bigger than 't_initial'
            The simulation's end time. The default value is 10.

        t_step : scalar, between 't_initial' and 't_final'
            The simulation's time step. The default value is 0.1.

        t_output : scalar, ND-array, or list
            When 't_output' is a scalar, it is considered as an output interval
            to store transient solutions. The default value is 1e+08. Initial,
            final and steady-state (if reached) fields are always stored. If
            't_output' > 't_final', no transient data is stored. If 't_output'
            is not a multiple of 't_step', 't_output' will be approximated.
            When 't_output' is a list or ND-array, transient solutions
            corresponding to this list or array will be stored.

        output_times : list
            List of output times. The values in the list must be multiples of
            the time step 't_step'.

        t_tolerance : scalar
            Transient solver tolerance. The simulation stops (before reaching
            't_final') when the residual falls below 't_tolerance'. The
            default value is 1e-06. The 'residual' measures the variation from
            one time-step to another in the value of the 'quantity' solved for.

        rxn_tolerance : scalar
            Tolerance to achieve within each time step. The solver passes to
            next time step when 'residual' falls below 'rxn_tolerance'. The
            default value is 1e-05.

        t_precision : integer
            The time precision (number of decimal places).

        t_scheme : string
            The time discretization scheme. Three options available: 'steady'
            to perform a steady-state simulation, and 'implicit' (fast, 1st
            order accurate) and 'cranknicolson' (slow, 2nd order accurate) both
            for transient simulations. The default value is 'implicit'.

        Notes
        -----
        More settings can be adjusted in the presence of a non-linear source
        term such as under-relaxation.
        See the 'ReactiveTransport' class documentation for details.
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
        if t_output is not None:
            self.settings['t_output'] = t_output
        if t_tolerance:
            self.settings['t_tolerance'] = t_tolerance
        if t_precision:
            self.settings['t_precision'] = t_precision
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
        # If ICs are not defined, assume zero
        try:
            self[self.settings['quantity']]
        except KeyError:
            self.set_IC(0)
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
        t_pre = self.settings['t_precision']
        s = self.settings['t_scheme']
        res_t = 1e+06  # Initialize the residual

        if type(to) in [float, int]:
            # Make sure 'tf' and 'to' are multiples of 'dt'
            tf = tf + (dt-(tf % dt))*((tf % dt) != 0)
            to = to + (dt-(to % dt))*((to % dt) != 0)
            self.settings['t_final'] = tf
            self.settings['t_output'] = to
            out = np.arange(t+to, tf, to)
        elif type(to) in [np.ndarray, list]:
            out = np.array(to)
        out = np.append(out, tf)
        out = np.unique(out)
        out = np.around(out, decimals=t_pre)

        if (s == 'steady'):  # If solver in steady mode, do one iteration
            logger.info('    Running in steady mode')
            x_old = self[self.settings['quantity']]
            self._t_run_reactive(x=x_old)
            x_new = self[self.settings['quantity']]

        else:  # Do time iterations
            # Export the initial field (t=t_initial)
            t_str = self._nbr_to_str(t)
            quant_init = self[self.settings['quantity']]
            self[self.settings['quantity']+'@'+t_str] = quant_init
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
                    if round(time, t_pre) in out:
                        t_str = self._nbr_to_str(time)
                        self[self.settings['quantity']+'@'+t_str] = x_new
                        logger.info('        Exporting time step: ' +
                                    str(time)+' s')
                    # Update A and b and apply BCs
                    self._t_update_A()
                    self._t_update_b()
                    self._apply_BCs()
                    self._A_t = (self._A).copy()
                    self._b_t = (self._b).copy()

                else:  # Stop time iterations if residual < t_tolerance
                    # Output steady state solution
                    t_str = self._nbr_to_str(time)
                    self[self.settings['quantity']+'@'+t_str] = x_new
                    logger.info('        Exporting time step: '+str(time)+' s')
                    break
            if (round(time, t_pre) == tf):
                logger.info('    Maximum time step reached: '+str(time)+' s')
            else:
                logger.info('    Transient solver converged after: ' +
                            str(time)+' s')

    def _t_run_reactive(self, x):
        """r
        Repeatedly updates transient 'A', 'b', and the solution guess within
        each time step according to the applied source term then calls '_solve'
        to solve the resulting system of linear equations. Stops when the
        residual falls below 'rxn_tolerance'.

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
        # Get the absolute values of non zero elements of A and b
        min_A = np.unique(np.abs(self.A.data))
        min_b = np.unique(np.abs(self.b[np.nonzero(self.b)]))
        # min of A & b after getting rid of the possible non conducting throats
        if min_A.size > 1:
            min_A = min_A[min_A != min_A.min()].min()
        else:
            min_A = 1
        if min_b.size > 1:
            min_b = min_b[min_b != min_b.min()].min()
        else:
            min_b = 1
        ref = min(min_A, min_b)  # Reference for residual's normalization
        for itr in range(int(self.settings['max_iter'])):
            self[self.settings['quantity']] = x
            self._A = (self._A_t).copy()
            self._b = (self._b_t).copy()
            self._apply_sources()
            # Compute the normalized residual
            res = np.linalg.norm(self.b-self.A*x)/ref
            if res >= self.settings['rxn_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                x_new = self._solve()
                # Relaxation
                x_new = relax*x_new + (1-relax)*self[self.settings['quantity']]
                self[self.settings['quantity']] = x_new
                x = x_new
            if (res < self.settings['rxn_tolerance'] or
                    self.settings['sources'] == []):
                x_new = x
                logger.info('Solution converged: ' + str(res))
                break
        return x_new

    def results(self, times=None, **kwargs):
        r"""
        Fetches the calculated quantity from the algorithm and returns it as
        an array.

        Parameters
        ----------
        times : scalar, ND-array, or list
            Time steps to be returned. The default value is None which results
            in returning all time steps. If times is a scalar, only the
            corresponding time step is returned. If times is an ND-array or a
            list, time steps in the provided array or list are returned.

        t_precision : integer
            The time precision (number of decimal places). Default value is 12.

        Notes
        -----
        The keyword steps is interpreted in the same way as times.
        """
        if 'steps' in kwargs.keys():
            times = kwargs['steps']
        t_pre = self.settings['t_precision']
        quantity = self.settings['quantity']
        q = [k for k in list(self.keys()) if quantity in k]
        if times is None:
            t = q
        elif type(times) in [np.ndarray, list, float, int]:
            out = np.array(times)
            out = np.unique(out)
            out = np.around(out, decimals=t_pre)
            t = []
            for i in out:
                j = self._nbr_to_str(i)
                t_str = [k for k in q if j == k.split('@')[-1]]
                t += (t_str)
            # Times stored by the transient algorithm
            strd_t = [float(st.split('@')[1]) for st in q if '@' in st]
            strd_t = np.array(strd_t)
            # Times requested but non stored by the algorithm
            missing_t = np.setdiff1d(np.around(out, decimals=t_pre),
                                     np.around(strd_t, decimals=t_pre))
            if missing_t.size != 0:
                logger.warning('Time(s) '+str(missing_t)+' not stored.')
        d = {k: self[k] for k in t}
        return d

    def _nbr_to_str(self, nbr, t_pre=None):
        r"""
        Converts a scalar into a string in scientific (exponential) notation
        without the decimal point.

        Parameters
        ----------
        nbr : scalar
            The number to be converted into a scalar.

        t_precision : integer
            The time precision (number of decimal places). Default value is 12.
        """
        if t_pre is None:
            t_pre = self.settings['t_precision']
        n = int(-dc(str(round(nbr, t_pre))).as_tuple().exponent *
                (round(nbr, t_pre) != int(nbr)))
        nbr_str = (str(int(round(nbr, t_pre)*10**n)) + ('e-'+str(n))*(n != 0))
        return nbr_str
