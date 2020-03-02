import numpy as np
from openpnm.algorithms import IonicTransport, TransientReactiveTransport


class TransientIonicTransport(IonicTransport, TransientReactiveTransport):
    r"""
    A subclass of GenericTransport to perform steady and transient simulations
    of pure diffusion, advection-diffusion and advection-diffusion with
    migration.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'gui': {'setup':        {'phase': None,
                                            'potential_field': '',
                                            'ions': [],
                                            'solver_tol': None,
                                            'solver_maxiter': None,
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_scheme': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, potential_field='', ions=[], i_tolerance=None,
              i_max_iter=None, t_initial=None, t_final=None, t_step=None,
              t_output=None, t_tolerance=None, t_precision=None, t_scheme='',
              **kwargs):
        if phase:
            self.settings['phase'] = phase.name
        if potential_field:
            self.settings['potential_field'] = potential_field
        if ions:
            self.settings['ions'] = ions
        if i_tolerance:
            self.settings['i_tolerance'] = i_tolerance
        if i_max_iter:
            self.settings['i_max_iter'] = i_max_iter
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
            try:
                alg[alg.settings['quantity']]
            except KeyError:
                try:
                    alg.set_IC(phase[alg.settings['quantity']])
                except KeyError:
                    alg.set_IC(0)

        # Save A matrix of the steady sys of eqs (WITHOUT BCs applied)
        for alg in algs:
            alg._build_A()
            alg._A_steady = (alg._A).copy()
        # Initialize A and b with BCs applied
        for e in e_alg:
            e._t_update_A()
            e._t_update_b()
            e._apply_BCs()
            e._A_t = (e._A).copy()
            e._b_t = (e._b).copy()
        # Init A&b with BCs for charge conservation eq, independent of t_scheme
        p_alg._apply_BCs()
        p_alg._A_t = (p_alg._A).copy()
        p_alg._b_t = (p_alg._b).copy()

        if t is None:
            t = self.settings['t_initial']
        # Create S1 & S1 for 1st Picard's iteration
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

        self._run_transient(t=t)

    def _run_transient(self, t):
        """r
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
        i_tol = self.settings['i_tolerance']
        t_pre = self.settings['t_precision']
        s = self.settings['t_scheme']
        # Initialize residuals & old/new fields for time marching
        t_res = {}
        t_old = {}
        t_new = {}
        for alg in algs:
            t_res[alg.name] = 1e+06
            t_old[alg.name] = None
            t_new[alg.name] = None

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

        # Source term for Poisson or charge conservation (electroneutrality) eq
        phys = p_alg.project.find_physics(phase=phase)
        p_alg._charge_conservation_eq_source_term(e_alg=e_alg)

        if (s == 'steady'):  # If solver in steady mode, do one iteration
            print('Running in steady mode')
            super().run()

        else:  # Do time iterations
            # Export the initial field (t=t_initial)
            t_str = self._nbr_to_str(t)
            for alg in algs:
                quant_init = alg[alg.settings['quantity']]
                alg[alg.settings['quantity']+'@'+t_str] = quant_init
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
                    i_res = {}
                    i_old = {}
                    i_new = {}
                    for alg in algs:
                        i_res[alg.name] = 1e+06
                        i_old[alg.name] = None
                        i_new[alg.name] = None

                    # Iterate (Gummel) until solutions converge
                    for itr in range(int(self.settings['i_max_iter'])):
                        i_r = [float(format(i, '.3g')) for i in i_res.values()]
                        i_r = str(i_r)[1:-1]
                        print('Gummel iter: '+str(itr+1)+', residuals: '+i_r)
                        i_convergence = max(i for i in i_res.values()) < i_tol
                        if not i_convergence:
                            # Ions
                            for e in e_alg:
                                i_old[e.name] = (
                                    e[e.settings['quantity']].copy())
                                e._t_run_reactive(x0=i_old[e.name])
                                i_new[e.name] = (
                                    e[e.settings['quantity']].copy())
                                # Residual
                                i_res[e.name] = np.sum(np.absolute(
                                    i_old[e.name]**2 - i_new[e.name]**2))
                                phase.update(e.results())

                            # Poisson eq
                            phys[0].regenerate_models()
                            i_old[p_alg.name] = (
                                p_alg[p_alg.settings['quantity']].copy())
                            p_alg._t_run_reactive(x0=i_old[p_alg.name])
                            i_new[p_alg.name] = (
                                p_alg[p_alg.settings['quantity']].copy())
                            # Residual
                            i_res[p_alg.name] = np.sum(np.absolute(
                                i_old[p_alg.name]**2 - i_new[p_alg.name]**2))
                            # Update phase and physics
                            phase.update(p_alg.results())
                            phys[0].regenerate_models()

                        elif i_convergence:
                            print('Solution for time step: ' + str(time)
                                  + ' s converged')
                            break

                    for alg in algs:  # Save new fields & compute t residuals
                        t_new[alg.name] = alg[alg.settings['quantity']].copy()
                        t_res[alg.name] = np.sum(
                            np.absolute(t_old[alg.name]**2 - t_new[alg.name]**2))

                    # Output transient solutions. Round time to ensure every
                    # value in outputs is exported.
                    if round(time, t_pre) in out:
                        t_str = self._nbr_to_str(time)
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
                        alg._A_steady = (alg._A).copy()

                    # Update A and b and apply BCs
                    for alg in algs:
                        alg._t_update_A()
                        alg._t_update_b()
                        alg._apply_BCs()
                        alg._A_t = (alg._A).copy()
                        alg._b_t = (alg._b).copy()

                else:  # Stop time iterations if residual < t_tolerance
                    # Output steady state solution
                    t_str = self._nbr_to_str(time)
                    print('\nExporting time step: '+str(time)+' s')
                    for alg in algs:
                        alg[alg.settings['quantity']+'@'+t_str] = (
                            t_new[alg.name])
                    break
            if (round(time, t_pre) == tf):
                print('\nMaximum time step reached: '+str(time)+' s')
            else:
                print('\nTransient solver converged after: '+str(time)+' s')
