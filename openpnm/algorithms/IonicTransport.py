import numpy as np
from openpnm.algorithms import ReactiveTransport


class IonicTransport(ReactiveTransport):
    r"""
    A subclass of GenericTransport to solve the charge conservation and
    Nernst-Planck equations.
    """
    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'potential_field': '',
                   'ions': [],
                   'i_tolerance': 1e-4,
                   'i_max_iter': 10}
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, potential_field='', ions=[], i_tolerance=None,
              i_max_iter=None, **kwargs):
        r"""
        """
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
        super().setup(**kwargs)

    def run(self, t=None):
        r"""
        """
        print('―'*80)
        print('Running IonicTransport')
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
                    alg[alg.settings['quantity']] = (
                        phase[alg.settings['quantity']])
                except KeyError:
                    alg[alg.settings['quantity']] = np.zeros(
                        shape=[alg.Np, ], dtype=float)

        # Source term for Poisson or charge conservation (electroneutrality) eq
        phys = p_alg.project.find_physics(phase=phase)
        p_alg._charge_conservation_eq_source_term(e_alg=e_alg)

        # Initialize residuals & old/new fields for Gummel iterats
        i_tol = self.settings['i_tolerance']
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
                    i_old[e.name] = (e[e.settings['quantity']].copy())
                    e._run_reactive(x0=i_old[e.name])
                    i_new[e.name] = (e[e.settings['quantity']].copy())
                    # Residual
                    i_res[e.name] = np.sum(np.absolute(
                        i_old[e.name]**2-i_new[e.name]**2))
                    phase.update(e.results())

                # Poisson eq
                phys[0].regenerate_models()
                i_old[p_alg.name] = p_alg[p_alg.settings['quantity']].copy()
                p_alg._run_reactive(x0=i_old[p_alg.name])
                i_new[p_alg.name] = p_alg[p_alg.settings['quantity']].copy()
                # Residual
                i_res[p_alg.name] = np.sum(np.absolute(
                    i_old[p_alg.name]**2 - i_new[p_alg.name]**2))
                # Update phase and physics
                phase.update(p_alg.results())
                phys[0].regenerate_models()

            if i_convergence:
                print('Solution converged')
                break
