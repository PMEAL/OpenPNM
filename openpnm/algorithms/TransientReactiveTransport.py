import numpy as np
import scipy.sparse as sprs
from decimal import Decimal as dc
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, GenericSettings, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sectionsf('TransientReactiveTransportSettings',
                      sections=['Parameters', 'Other Parameters'])
@docstr.dedent
class TransientReactiveTransportSettings(GenericSettings):
    r"""

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s

    quantity : (str)
        The name of the physical quantity to be calculated
    conductance : (str)
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------
    t_initial : scalar
        The simulation's start time, which must be smaller than
        't_final'. The default value is 0.
    t_final : scalar
        The simulation's end time, which must be bigger than 't_initial'.
        The default value is 10.
    t_step : scalar
        The simulation's time step, which must be smaller than 't_final' less
        and 't_initial'. The default value is 0.1.
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
    t_precision : integer
        The time precision (number of decimal places).
    t_scheme : string
        The time discretization scheme. Three options available: 'steady'
        to perform a steady-state simulation, and 'implicit' (fast, 1st
        order accurate) and 'cranknicolson' (slow, 2nd order accurate) both
        for transient simulations. The default value is 'implicit'.

    ----

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """

    phase = None
    t_initial = 0
    t_final = 10
    t_step = 0.1
    t_output = 1e+08
    t_tolerance = 1e-06
    t_precision = 12
    t_scheme = 'implicit'


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
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(TransientReactiveTransportSettings)
        self.settings.update(settings)
        # Initialize the steady sys of eqs A matrix
        self._A_steady = None
        if phase is not None:
            self.setup(phase=phase)
        # Initialize the initial condition
        self["pore.ic"] = np.nan

    def setup(self, phase=None, quantity='', conductance='',
              t_initial=None, t_final=None, t_step=None, t_output=None,
              t_tolerance=None, t_precision=None, t_scheme='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings

        Parameters
        ----------

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
        if t_initial is not None:
            self.settings['t_initial'] = t_initial
        if t_final is not None:
            self.settings['t_final'] = t_final
        if t_step is not None:
            self.settings['t_step'] = t_step
        if t_output is not None:
            self.settings['t_output'] = t_output
        if t_tolerance is not None:
            self.settings['t_tolerance'] = t_tolerance
        if t_precision is not None:
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
        quantity = self.settings['quantity']
        values = np.array(values) * 1.0
        if values.size > 1 and values.size != self.Np:
            raise Exception('The number of initial values must be either 1 or Np')
        self['pore.ic'] = values
        if not quantity:
            raise Exception('"quantity" has not been defined on this algorithm')
        self[quantity] = values

    def _get_f1_f2_f3(self):
        r"""
        Helper method: returns f1, f2, and f3 for _t_update_A and _t_update_b methods.
        """
        s = self.settings['t_scheme']
        if s == 'implicit':
            f1, f2, f3 = 1, 1, 0
        elif s == 'cranknicolson':
            f1, f2, f3 = 0.5, 1, 0
        elif s == 'steady':
            f1, f2, f3 = 1, 0, 1
        else:
            raise Exception(f'Unsupported t_scheme: "{s}"')
        return f1, f2, f3

    def _t_update_A(self):
        r"""
        A method to update 'A' matrix at each time step according to 't_scheme'
        """
        network = self.project.network
        Vi = network['pore.volume']
        dt = self.settings['t_step']
        f1, f2, _ = self._get_f1_f2_f3()
        # Compute A (operations involve conversion to 'csr')
        A = ((f2/dt) * sprs.coo_matrix.multiply(
            sprs.coo_matrix(np.reshape(Vi, (self.Np, 1)), shape=(self.Np, 1)),
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
        quantity = self.settings['quantity']
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        Vi = network['pore.volume']
        dt = self.settings['t_step']
        f1, f2, f3 = self._get_f1_f2_f3()
        x_old = self[quantity]
        b = (f2 * (1-f1) * (-self._A_steady) * x_old
             + f2 * (Vi/dt) * x_old
             + f3 * np.zeros(shape=(self.Np,), dtype=float))
        self._update_iterative_props()
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Update b
            b[Ps] = b[Ps] - f2 * (1-f1) * (phase[item + '.' + 'rate'][Ps])
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
        logger.info('―' * 80)
        logger.info('Running TransientTransport')
        # If ICs are not defined, assume zero
        if not np.isfinite(self["pore.ic"]).all():
            self.set_IC(0)
        # Make sure _A is None to force _build_A, otherwise _A_steady might be wrong
        self._A = None
        # Save A matrix of the steady state problem (without BCs applied)
        self._A_steady = self.A.copy()
        # Initialize A and b with BCs applied
        self._t_update_A()
        self._t_update_b()
        self._apply_BCs()
        # Save copies of A and b to be used in _t_run_reactive()
        self._A_t = self._A.copy()
        self._b_t = self._b.copy()
        if t is None:
            t = self.settings['t_initial']
        # Create S1 & S2 for 1st Picard's iteration
        self._update_iterative_props()
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
        quantity = self.settings['quantity']
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

        # If solver in steady mode, do one iteration
        if s == 'steady':
            logger.info('    Running in steady mode')
            self._t_run_reactive()

        # Time marching step
        else:
            # Export the initial field (t=t_initial)
            t_str = self._nbr_to_str(t)
            quant_init = self["pore.ic"]
            self[quantity + '@' + t_str] = quant_init
            self[quantity] = quant_init

            for time in np.arange(t+dt, tf+dt, dt):
                if res_t >= tol:  # Check if the steady state is reached
                    logger.info(f'    Current time step: {time} s')
                    # Update A and b and apply BCs
                    self._t_update_A()
                    self._t_update_b()
                    self._apply_BCs()
                    # Save copies of A and b to be used in _t_run_reactive()
                    self._A_t = self._A.copy()
                    self._b_t = self._b.copy()
                    x_old = self[quantity]
                    self._t_run_reactive(x0=x_old)
                    x_new = self[quantity]
                    # Compute the residual
                    res_t = np.sum(np.absolute(x_old**2 - x_new**2))
                    logger.info(f'        Residual: {res_t}')
                    # Output transient solutions. Round time to ensure every
                    # value in outputs is exported.
                    if round(time, t_pre) in out:
                        t_str = self._nbr_to_str(time)
                        self[quantity + '@' + t_str] = x_new
                        logger.info(f'        Exporting time step: {time} s')
                else:
                    # Output steady state solution
                    t_str = self._nbr_to_str(time)
                    self[quantity + '@' + t_str] = x_new
                    logger.info(f'        Exporting time step: {time} s')
                    break

            if round(time, t_pre) == tf:
                logger.info(f'    Maximum time step reached: {time} s')
            else:
                logger.info(f'    Transient solver converged after: {time} s')

    def _t_run_reactive(self, x0=None):
        """r
        Repeatedly updates transient 'A', 'b', and the solution guess within
        each time step according to the applied source term then calls '_solve'
        to solve the resulting system of linear equations. Stops when the
        residual falls below 'solver_tol'.

        Parameters
        ----------
        x0 : ND-array
            Initial guess of unknown variable

        Returns
        -------
        x_new : ND-array
            Solution array.

        Notes
        -----
        Description of 'relaxation_quantity' and 'nlin_max_iter' settings can be
        found in the parent class 'ReactiveTransport' documentation.

        """
        quantity = self.settings['quantity']
        w = self.settings['relaxation_quantity']
        max_it = int(self.settings['nlin_max_iter'])
        x = np.zeros(self.Np, dtype=float) if x0 is None else x0.copy()

        # Write initial guess to algorithm for _update_iterative_props to work
        self[quantity] = x

        for itr in range(max_it):
            # Update iterative properties on phase and physics
            self._update_iterative_props()
            # Build A and b, apply source terms and correct according to scheme
            self._A = self._A_t.copy()
            self._b = self._b_t.copy()
            self._apply_sources()
            self._correct_apply_sources()
            # Compute the residual
            res = self._get_residual()
            if itr >= 1 and self._is_converged():
                logger.info(f'Solution converged: {res:.4e}')
                return x
            logger.info(f'Tolerance not met: {res:.4e}')
            # Solve, use relaxation, and update solution on algorithm obj
            self[quantity] = x = self._solve(x0=x) * w + x * (1 - w)

        # Check solution convergence after max_it iterations
        if not self._is_converged():
            raise Exception(f"Not converged after {max_it} iterations.")

    def results(self, times=None, **kwargs):
        r"""
        Fetches the calculated quantity from the algorithm and returns it as
        an array.

        Parameters
        ----------
        times : scalar, ND-array, list of scalars, None, or string
            Time steps to be returned. The default value is None which results
            in returning all time steps. If times is a scalar, only the
            corresponding time step is returned. If times is an ND-array or a
            list of scalars, time steps in the provided array or list are
            returned. If times is 'final' or 'actual', the current value of the
            quantity is returned.

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
        elif times in ['final', 'actual']:
            t = [quantity]
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
        n = int(-dc(str(round(nbr, t_pre))).as_tuple().exponent
                * (round(nbr, t_pre) != int(nbr)))
        nbr_str = (str(int(round(nbr, t_pre)*10**n)) + ('e-'+str(n))*(n != 0))
        return nbr_str

    def _correct_apply_sources(self):
        """r
        Update 'A' and 'b' correcting the already applied source terms to
        specified pores

        Notes
        -----
        Correction (built for transient simulations) depends on the time scheme
        """
        if self.settings['t_scheme'] == 'cranknicolson':
            f1 = 0.5
        else:
            f1 = 1.0
        phase = self.project.phases()[self.settings['phase']]
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # get already added relaxed source term
            S1, S2 = [phase[item + '.' + x][Ps] for x in ['S1', 'S2']]
            # correct S1 and S2 in A and b as a function of t_scheme
            datadiag = self._A.diagonal().copy()
            datadiag[Ps] = datadiag[Ps] - S1 + f1*S1
            self._A.setdiag(datadiag)
            self._b[Ps] = self._b[Ps] + S2 - f1*S2
