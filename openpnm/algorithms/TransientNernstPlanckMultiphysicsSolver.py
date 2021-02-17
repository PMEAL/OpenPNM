import numpy as np
from openpnm.algorithms import NernstPlanckMultiphysicsSolver
from openpnm.utils import logging, Docorator, GenericSettings, nbr_to_str
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='TransientNernstPlanckMultiphysicsSolverSettings',
                     sections=['Parameters'])
@docstr.dedent
class TransientNernstPlanckMultiphysicsSolverSettings(GenericSettings):
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

    %(NernstPlanckMultiphysicsSolverSettings.parameters)s

    """
    t_initial = 0
    t_final = 10
    t_step = 0.1
    t_output = 1e+08
    t_tolerance = 1e-06
    t_precision = 12
    t_scheme = 'implicit'


class TransientNernstPlanckMultiphysicsSolver(NernstPlanckMultiphysicsSolver):
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
        c = TransientNernstPlanckMultiphysicsSolverSettings()
        self.settings._update_settings_and_docs(c)
        self.settings.update(settings)

    @docstr.dedent
    def setup(self, t_initial=None, t_final=None, t_step=None,
              t_output=None, t_tolerance=None, t_precision=None,
              t_scheme='implicit', **kwargs):
        r"""

        Parameters
        ----------
        %(TransientNernstPlanckMultiphysicsSolverSettings.parameters)s

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
        # Phase, potential and ions algorithms
        phase = self.project.phases()[self.settings['phase']]
        p_alg = self.project.algorithms()[self.settings['potential_field']]
        e_alg = [self.project.algorithms()[self.settings['ions'][i]] for i in
                 range(len(self.settings['ions']))]
        algs = e_alg.copy()
        algs.insert(0, p_alg)
        # Define initial conditions (if not defined by the user)
        for alg in algs:
            alg.settings.update({'cache_A': False, 'cache_b': False})
            try:
                alg[alg.settings['quantity']]
            except KeyError:
                try:
                    alg.set_IC(phase[alg.settings['quantity']])
                except KeyError:
                    alg.set_IC(0)

        for e in e_alg:
            # Save A matrix of the steady sys of eqs (WITHOUT BCs applied)
            e._build_A()
            e._A_steady = e._A.copy()
            # Initialize A and b with BCs applied
            e._t_update_A()
            e._t_update_b()
            e._apply_BCs()
            e._A_t = e._A.copy()
            e._b_t = e._b.copy()
        # Init A&b with BCs for charge conservation eq, independent of t_scheme
        p_alg._build_A()
        p_alg._apply_BCs()

        if t is None:
            t = self.settings['t_initial']
        # Create S1 & S1 for 1st Picard's iteration
        for alg in algs:
            alg._update_iterative_props()

        # Setup algorithms transient settings
        for e in e_alg:
            e.setup(t_initial=self.settings['t_initial'],
                    t_final=self.settings['t_final'],
                    t_step=self.settings['t_step'],
                    t_output=self.settings['t_output'],
                    t_tolerance=self.settings['t_tolerance'],
                    t_precision=self.settings['t_precision'],
                    t_scheme=self.settings['t_scheme'])

        self._run_transient(t=t)

    def _run_transient(self, t):
        r"""

        """
        # Phase, potential and ions algorithms
        phase = self.project.phases()[self.settings['phase']]
        p_alg = self.project.algorithms()[self.settings['potential_field']]
        e_alg = [self.project.algorithms()[self.settings['ions'][i]] for i in
                 range(len(self.settings['ions']))]
        algs = e_alg.copy()
        algs.insert(0, p_alg)

        tf = self.settings['t_final']
        dt = self.settings['t_step']
        to = self.settings['t_output']
        t_tol = self.settings['t_tolerance']
        t_pre = self.settings['t_precision']
        s = self.settings['t_scheme']
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

        # Source term for Poisson or charge conservation (electroneutrality) eq
        phys = p_alg.project.find_physics(phase=phase)
        p_alg._charge_conservation_eq_source_term(e_alg=e_alg)

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
                            # Ions
                            for e in e_alg:
                                g_old[e.name] = (
                                    e[e.settings['quantity']].copy())

                                e._t_run_reactive(x0=g_old[e.name])
                                g_new[e.name] = (
                                    e[e.settings['quantity']].copy()
                                )
                                # Residual
                                g_res[e.name] = np.sum(np.absolute(
                                    g_old[e.name]**2 - g_new[e.name]**2))
                                phase.update(e.results())

                            # Charge conservation eq
                            for obj in phys:
                                obj.regenerate_models()
                            g_old[p_alg.name] = (
                                p_alg[p_alg.settings['quantity']].copy())
                            p_alg._run_reactive(x0=g_old[p_alg.name])
                            g_new[p_alg.name] = (
                                p_alg[p_alg.settings['quantity']].copy()
                            )
                            # Residual
                            g_res[p_alg.name] = np.sum(np.absolute(
                                g_old[p_alg.name]**2 - g_new[p_alg.name]**2))
                            # Update phase and physics
                            phase.update(p_alg.results())
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
                    for e in e_alg:
                        # Update conductance first
                        physics = e.project.find_physics(phase=phase)
                        for ph in physics:
                            ph.regenerate_models()
                        # Update A matrix
                        e._build_A()
                        e._A_steady = e._A.copy()

                    # Update A and b and apply BCs
                    for e in e_alg:
                        e._t_update_A()
                        e._t_update_b()
                        e._apply_BCs()
                        e._A_t = e._A.copy()
                        e._b_t = e._b.copy()

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
