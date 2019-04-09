import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
from openpnm.models.physics import generic_source_term as gst
logger = logging.getLogger(__name__)


class PoissonNernstPlanck(ReactiveTransport):
    r"""
    A subclass of GenericTransport to solve the Poisson Nernst-Planck equations
    """
    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'potential_field': None,
                   'electrolytes': [],
                   'charge_conservation': 'electroneutrality',
                   'tolerance': 1e-4,
                   'max_iter': 10}
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, potential_field=None, electrolytes=[],
              charge_conservation=None, tolerance=None,
              max_iter=None, **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if potential_field:
            self.settings['potential_field'] = potential_field
        if electrolytes:
            self.settings['electrolytes'] = electrolytes
        if charge_conservation:
            self.settings['charge_conservation'] = charge_conservation
        if tolerance:
            self.settings['tolerance'] = tolerance
        if max_iter:
            self.settings['max_iter'] = max_iter
        super().setup(**kwargs)

    def run(self):
        # Phase, potential and electrolytes algorithms
        phase = self.project.phases()[self.settings['phase']]
        p_alg = self.settings['potential_field']
        e_alg = self.settings['electrolytes']

        # Define initial conditions (if not defined by the user)
        try:
            p_alg[p_alg.settings['quantity']]
        except KeyError:
            p_alg[p_alg.settings['quantity']] = np.zeros(shape=[p_alg.Np, ],
                                                         dtype=float)
        for e in e_alg:
            try:
                e[e.settings['quantity']]
            except KeyError:
                e[e.settings['quantity']] = np.zeros(shape=[e.Np, ],
                                                     dtype=float)

        # Define source term for Poisson equation
        Ps = (p_alg['pore.all'] * np.isnan(p_alg['pore.bc_value']) *
              np.isnan(p_alg['pore.bc_rate']))
        mod = gst.charge_conservation
        phys = p_alg.project.find_physics(phase=phase)
        phys[0].add_model(propname='pore.charge_conservation', model=mod,
                          phase=phase, p_alg=p_alg, e_alg=e_alg,
                          assumption=self.settings['charge_conservation'])
        p_alg.set_source(propname='pore.charge_conservation', pores=Ps)

        # Define tolerance and initialize residuals
        tol = self.settings['tolerance']
        res = {}
        res['potential'] = 1e+06
        for e in e_alg:
            res[e.name] = 1e+06

        # Iterate until solutions converge
        for itr in range(int(self.settings['max_iter'])):
            r = str([float(format(i, '.3g')) for i in res.values()])[1:-1]
            logger.info('Iter: ' + str(itr) + ', Residuals: ' + r)
            print('Iter: ' + str(itr) + ', Residuals: ' + r)
            convergence = max(i for i in res.values()) < tol
            if not convergence:
                # Poisson eq
                phys[0].regenerate_models()
                phi_old = p_alg[p_alg.settings['quantity']].copy()
                p_alg._run_reactive(x=phi_old)
                phi_new = p_alg[p_alg.settings['quantity']].copy()
                # Residual
                res['potential'] = np.sum(np.absolute(phi_old**2 - phi_new**2))
                # Update phase and physics
                phase.update(p_alg.results())
                phys[0].regenerate_models()

                # Electrolytes
                for e in e_alg:
                    c_old = e[e.settings['quantity']].copy()
                    e._run_reactive(x=c_old)
                    c_new = e[e.settings['quantity']].copy()
                    # Residual
                    res[e.name] = np.sum(np.absolute(c_old**2 - c_new**2))
                    phase.update(e.results())

            if convergence:
                logger.info('Solution converged: ' + str(res))
                break
