import numpy as np
from openpnm.algorithms import MultiphysicsSolver
from openpnm.utils import logging, Docorator, GenericSettings, nbr_to_str
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sectionsf('MultiphysicsSolverSettings',
                      sections=['Parameters'])
@docstr.dedent
class TransientMultiphysicsSolverSettings(GenericSettings):
    r"""
    The Parameters section below describes the settings pertaining to the
    running of all transient classes which this algorithm orchestrates.

    Parameters
    ----------
    %(TransientReactiveTransportSettings.other_parameters)s

    Other Parameters
    ----------------
    **The following parameters pertain to the steady-state version of this
    class**

    %(MultiphysicsSolverSettings.parameters)s

    """
    t_initial = 0
    t_final = 10
    t_step = 0.1
    t_output = 1e+08
    t_tolerance = 1e-06
    t_precision = 12
    t_scheme = 'implicit'


class TransientMultiphysicsSolver(MultiphysicsSolver):
    r"""
    A multiphysics solver to solve the Nernst-Planck and Ionic Conduction
    system *transiently*.

    Warnings
    --------
    This is not a true OpenPNM algorithm. This solver wraps the provided
    Nernst-Planck and ionic conduction algorithms and solves the associated
    system of equations.

    """
    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        c = TransientMultiphysicsSolverSettings()
        self.settings._update_settings_and_docs(c)
        self.settings.update(settings)

    @docstr.dedent
    def setup(self, t_initial=None, t_final=None, t_step=None,
              t_output=None, t_tolerance=None, t_precision=None,
              t_scheme='implicit', **kwargs):
        r"""

        Parameters
        ----------
        %(TransientMultiphysicsSolverSettings.parameters)s

        """
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
        self.settings.update(**kwargs)

    def run(self, t=None):
        r"""

        """
        print('―'*80)
        print('Running TransientIonicTransport')
        # Retrieve algorithm objects
        algs = self._get_algorithms()
        # Force cache_A and cache_B to False
        self._set_cache_to_false(algorithms=algs)
        # Initialize transient A and b matrices
        for alg in algs:
            # Save A matrix of the steady sys of eqs (WITHOUT BCs applied)
            alg._build_A()
            alg._A_steady = alg._A.copy()
            # Initialize A and b with BCs applied.
            alg._t_update_A()
            alg._t_update_b()
            alg._apply_BCs()
            alg._A_t = alg._A.copy()
            alg._b_t = alg._b.copy()
        # get initial time
        if t is None:
            t = self.settings['t_initial']
        # Create S1 & S2 for 1st Picard's iteration
        for alg in algs:
            alg._update_iterative_props()
        # Setup algorithms transient settings
        for alg in algs:
            alg.setup(t_initial=self.settings['t_initial'],
                      t_final=self.settings['t_final'],
                      t_step=self.settings['t_step'],
                      t_output=self.settings['t_output'],
                      t_tolerance=self.settings['t_tolerance'],
                      t_precision=self.settings['t_precision'],
                      t_scheme=self.settings['t_scheme'])
        # run transient multiphysics solver
        self._run_transient(t=t)

    def _run_transient(self, t):
        r"""

        """
        # Retrieve phase, physics, and algorithm objects
        phase = self.project.phases()[self.settings['phase']]
        phys = self.project.find_physics(phase=phase)
        algs = self._get_algorithms()
        # Retrieve transient settings from solver
        tf = self.settings['t_final']
        dt = self.settings['t_step']
        to = self.settings['t_output']
        t_tol = self.settings['t_tolerance']
        t_pre = self.settings['t_precision']
        s = self.settings['t_scheme']
        # retrieve gummel settings from solver
        g_tol = self.settings['g_tol']
        g_max_iter = int(self.settings['g_max_iter'])
        # Initialize residuals & old/new fields for time marching
        t_res = {}
        t_old = {}
        t_new = {}
        for alg in algs:
            t_res[alg.name] = 1e+06
            t_old[alg.name] = None
            t_new[alg.name] = None

        if isinstance(to, (float, int)):
            # Make sure 'tf' and 'to' are multiples of 'dt'
            tf = tf + (dt-(tf % dt))*((tf % dt) != 0)
            to = to + (dt-(to % dt))*((to % dt) != 0)
            self.settings['t_final'] = tf
            self.settings['t_output'] = to
            out = np.arange(t+to, tf, to)
        elif isinstance(to, (np.ndarray, list)):
            out = np.array(to)
        out = np.append(out, tf)
        out = np.unique(out)
        out = np.around(out, decimals=t_pre)

        if s == 'steady':  # If solver in steady mode, do one iteration
            print('Running in steady mode')
            super().run()

        else:  # Do time iterations
            # Export the initial field (t=t_initial)
            t_str = nbr_to_str(nbr=t, t_precision=self.settings['t_precision'])
            for alg in algs:
                quant_init = alg[alg.settings['quantity']]
                alg[alg.settings['quantity']+'@'+t_str] = quant_init
            time = t + dt
            for time in np.arange(t+dt, tf+dt, dt):
                t_r = [float(format(i, '.3g')) for i in t_res.values()]
                t_r = str(t_r)[1:-1]
                print('\n'+'Current time step: '+str(time)+' s')
                print('Algorithms: '+', '.join(t_res.keys()))
                print('Time residuals: '+t_r)
                t_convergence = max(i for i in t_res.values()) < t_tol
                if not t_convergence:  # Check if the steady state is reached
                    for alg in algs:  # Save the current fields
                        t_old[alg.name] = alg[alg.settings['quantity']].copy()

                    # Initialize residuals & old/new fields for Gummel iterats
                    g_res = {}
                    g_old = {}
                    g_new = {}
                    for alg in algs:
                        g_res[alg.name] = 1e+03
                        g_old[alg.name] = None
                        g_new[alg.name] = None

                    # Iterate (Gummel) until solutions converge
                    for itr in range(g_max_iter):
                        g_r = [float(format(i, '.3g')) for i in g_res.values()]
                        g_r = str(g_r)[1:-1]
                        print('Start Gummel iter: ' + str(itr+1)
                              + ', residuals: ' + g_r)
                        g_convergence = max(i for i in g_res.values()) < g_tol
                        if not g_convergence:
                            # Solve each algorithm consecutively
                            for alg in algs:
                                g_old[alg.name] = (alg[alg.settings['quantity']].copy())
                                alg._t_run_reactive(x0=g_old[alg.name])
                                g_new[alg.name] = (alg[alg.settings['quantity']].copy())
                                # Residual
                                g_res[alg.name] = np.sum(np.absolute(
                                    g_old[alg.name]**2-g_new[alg.name]**2))
                                # update phase and regenerate physics
                                phase.update(alg.results())
                                for obj in phys:
                                    obj.regenerate_models()

                        elif g_convergence:
                            print('Solution for time step: ' + str(time)
                                  + ' s converged')
                            break

                    for alg in algs:  # Save new fields & compute t residuals
                        t_new[alg.name] = alg[alg.settings['quantity']].copy()
                        t_res[alg.name] = np.sum(
                            np.absolute(t_old[alg.name]**2
                                        - t_new[alg.name]**2))

                    # Output transient solutions. Round time to ensure every
                    # value in outputs is exported.
                    if round(time, t_pre) in out:
                        t_str = nbr_to_str(nbr=time,
                                           t_precision=self.settings['t_precision'])
                        print('\nExporting time step: ' + str(time) + ' s')
                        for alg in algs:
                            alg[alg.settings['quantity']+'@'+t_str] = (
                                t_new[alg.name])
                        
                    # Update A matrix of the steady sys of eqs (WITHOUT BCs)
                    for alg in algs:
                        # Update conductance first
                        physics = alg.project.find_physics(phase=phase)
                        for ph in physics:
                            ph.regenerate_models()
                        # Update A matrix
                        alg._build_A()
                        alg._A_steady = alg._A.copy()
                        # Update A and b and apply BCs
                        alg._t_update_A()
                        alg._t_update_b()
                        alg._apply_BCs()
                        alg._A_t = alg._A.copy()
                        alg._b_t = alg._b.copy() 

                    # update _old properties
                    for obj in phys:
                        for alg in algs:
                            old_prop = alg.settings['old_prop']
                            if old_prop is not None:
                                obj[old_prop + '_old'] = obj[old_prop].copy()
                            
                else:  # Stop time iterations if residual < t_tolerance
                    # Output steady state solution
                    t_str = nbr_to_str(nbr=time,
                                       t_precision=self.settings['t_precision'])
                    print('\nExporting time step: '+str(time)+' s')
                    for alg in algs:
                        alg[alg.settings['quantity']+'@'+t_str] = (
                            t_new[alg.name])
                    break
            if round(time, t_pre) == tf:
                print('\nMaximum time step reached: '+str(time)+' s')
            else:
                print('\nTransient solver converged after: '+str(time)+' s')
