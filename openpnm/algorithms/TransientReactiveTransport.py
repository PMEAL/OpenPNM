import numpy as np
import scipy.sparse as sprs
from openpnm.algorithms import ReactiveTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class TransientReactiveTransport(ReactiveTransport):
    r"""
    A subclass of GenericTransport to perform transient and steady simulations.
    """

    def __init__(self, settings={}, **kwargs):
        self.settings.update({'t_initial': 0,
                              't_final': 5e+04,
                              't_step': 0.1,
                              't_output': 1e+08,
                              't_tolerance': 1e-06,
                              'r_tolerance': 1e-04,
                              't_scheme': 'implicit'})
        super().__init__(**kwargs)
        self.settings.update(settings)
        self._A_steady = None  # Initialize the steady sys of eqs A matrix

    def setup(self, phase=None, t_initial='', t_final='', t_step='',
              t_output='', t_tolerance='', t_scheme='', **kwargs):
        if phase:
            self.settings['phase'] = phase.name
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
        """
        self[self.settings['quantity']] = values

    def _t_update_A(self):
        r"""
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
        """
        print('―'*80)
        print('Running TransientTransport')
        # If solver used in steady mode, no need to add ICs
        if (self.settings['t_scheme'] == 'steady'):
            self[self.settings['quantity']] = 0
        # Create a scratch b from IC
        self._b = (self[self.settings['quantity']]).copy()
        self._apply_BCs()
        # Save A matrix (with BCs applied) of the steady sys of eqs
        self._A_steady = (self._A).copy()
        # Save the initial field with the boundary conditions applied
        self[self.settings['quantity']] = (self._b).copy()
        # Override A and b according to t_scheme and apply BCs
        self._t_update_A()
        self._t_update_b()
        self._apply_BCs()
        self._A_t = (self._A).copy()
        self._b_t = (self._b).copy()
        if t is None:
            t = self.settings['t_initial']
        self._run_transient(t=t)

    def _run_transient(self, t):
        tf = self.settings['t_final']
        dt = self.settings['t_step']
        to = self.settings['t_output']
        tol = self.settings['t_tolerance']
        s = self.settings['t_scheme']
        res_t = 1  # Initialize the residual

        # Make sure 'tf' and 'to' are multiples of 'dt'
        tf = tf + (dt-(tf % dt))*((tf % dt) != 0)
        to = to + (dt-(to % dt))*((to % dt) != 0)
        self.settings['t_final'] = tf
        self.settings['t_output'] = to
        outputs = np.append(np.arange(t+to, tf, to), tf)

        if (s == 'steady'):  # If solver in steady mode, do one iteration
            print('    Running in steady mode')
            x_old = self[self.settings['quantity']]
            self._t_run_reactive(x=x_old)
            x_new = self[self.settings['quantity']]

        else:  # Do time iterations
            # Export the initial field (t=t_initial)
            vals = self[self.settings['quantity']]
            self[self.settings['quantity']+'_initial'] = vals
            for time in np.arange(t+dt, tf+dt, dt):
                if (res_t >= tol):  # Check if the steady state is reached
                    print('    Current time step: '+str(time)+' s')
                    x_old = self[self.settings['quantity']]

                    self._t_run_reactive(x=x_old)
                    x_new = self[self.settings['quantity']]

                    # Compute the residual
                    res_t = np.sum(np.absolute(x_old**2 - x_new**2))
                    print('        Residual: '+str(res_t))
                    # Output transient solutions. Round time to ensure every
                    # value in outputs is exported.
                    if round(time, 12) in outputs:
                        ind = np.where(outputs == round(time, 12))[0][0]
                        self[self.settings['quantity']+'_'+str(ind)] = x_new
                        print('        Exporting time step: '+str(time)+' s')
                    # Update b and apply BCs
                    self._t_update_b()
                    self._apply_BCs()
                    self._b_t = (self._b).copy()
                else:  # Stop time iterations if residual < t_tolerance
                    self[self.settings['quantity'] + '_steady'] = x_new
                    print('        Exporting time step: '+str(time)+' s')
                    break
            if (round(time, 12) == tf):
                print('    Maximum time step reached: '+str(time)+' s')
            else:
                print('    Transient solver converged after: '+str(time)+' s')

    def _t_run_reactive(self, x):
        if self.settings['quantity'] not in self.keys():
            self[self.settings['quantity']] = 0
        self._A = (self._A_t).copy()
        self._b = (self._b_t).copy()
        self._apply_BCs()
        self._t_apply_sources(s=self.settings['t_scheme'])
        if x is None:
            x = np.zeros(shape=[self.Np, ], dtype=float)
        x_new = self._solve()
        self[self.settings['quantity']] = x_new
        res_r = np.sum(np.absolute(x**2 - x_new**2))
        if res_r < self.settings['r_tolerance']:
            print('            Solution converged: ' + str(res_r))
            return x_new
        else:
            print('            Tolerance not met: ' + str(res_r))
            self._t_run_reactive(x=x_new)

    def _t_apply_sources(self, s=None):
        f1 = 1
        if (s == 'cranknicolson'):
            f1 = 0.5
        phase = self.project.phases()[self.settings['phase']]
        self._update_physics()
        for item in self.settings['sources']:
            Ps = self.pores(item)
            # Add S1 to diagonal of A
            # TODO: We need this to NOT overwrite the A and b, but create
            # copy, otherwise we have to regenerate A and b on each loop
            datadiag = self.A.diagonal()
            datadiag[Ps] = datadiag[Ps] + f1*phase[item+'.'+'S1'][Ps]
            self.A.setdiag(datadiag)
            # Add S2 to b
            self.b[Ps] = self.b[Ps] - f1*phase[item+'.'+'S2'][Ps]
