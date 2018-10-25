import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class PoissonNernstPlanck(ReactiveTransport):
    r"""
    A subclass of GenericTransport to solve the Poisson Nernst-Planck equations

    """
    def __init__(self, settings={}, phase=None, poisson_alg=None, **kwargs):
        def_set = {'phase': None,
                   'poisson_alg': poisson_alg,
                   'hydraulic_conductance': 'throat.hydraulic_conductance',
                   'pressure': 'pore.pressure',
                   'solvent': {'permittivity': None},
                   'electrolytes': {'electrolyteA': {'algorithm': None,
                                                     'valence': None,
                                                     'diffusion_coef': None},
                                    'electrolyteB': {'algorithm': None,
                                                     'valence': None,
                                                     'diffusion_coef': None}},
                   'pnp_tolerance': 1e-06}
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def _build_A(self, electrolyte=None, force=False):
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        conns = network['throat.conns']

        P = phase[self.settings['pressure']]
        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase['throat.diffusive_conductance'+'_'+str(electrolyte)]
        gd = np.tile(gd, 2)

        Qij = -gh*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)

        Peij = Qij/gd
        Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
        Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10
        Qij = Peij*gd

        p_alg = self.settings['poisson_alg']
        phi = p_alg[p_alg.settings['quantity']]
        S = network['throat.area']
        L = network['throat.length']
        z = self.settings['electrolytes'][electrolyte]['valence']
        D = self.settings['electrolytes'][electrolyte]['diffusion_coef']
        F = 96485.3329
        R = 8.3145
        T = 298

        grad_phi = -np.diff(phi[conns], axis=1).squeeze() / L
        mig = ((z*F*D*S)/(R*T)) * grad_phi
        mig = np.append(mig, -mig)

        flow_mig = Qij-mig

        if force:
            self._pure_A = None
        if self._pure_A is None:
            w = gd + np.maximum(0, -flow_mig)
            A = network.create_adjacency_matrix(weights=w)
            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()
        return self.A

    def _update_b_poisson(self):
        p_alg = self.settings['poisson_alg']
        eA_alg = self.settings['electrolytes']['electrolyteA']['algorithm']
        eB_alg = self.settings['electrolytes']['electrolyteB']['algorithm']
        F = 96485.3329
        epsilon = self.settings['solvent']['permittivity']
        zA = self.settings['electrolytes']['electrolyteA']['valence']
        zB = self.settings['electrolytes']['electrolyteB']['valence']
        cA = eA_alg[eA_alg.settings['quantity']]
        cB = eB_alg[eB_alg.settings['quantity']]
        p_alg.b = (-F/epsilon) * (zA*cA + zB*cB)

    def run(self):
        p_alg = self.settings['poisson_alg']
        eA_alg = self.settings['electrolytes']['electrolyteA']['algorithm']
        eB_alg = self.settings['electrolytes']['electrolyteB']['algorithm']
        p_alg[p_alg.settings['quantity']] = 0
        eA_alg[eA_alg.settings['quantity']] = 0
        eB_alg[eB_alg.settings['quantity']] = 0

        # Initialize the residuals
        resP = 1e+06
        resA = 1e+06
        resB = 1e+06

        for itr in range(int(self.settings['max_iter'])):
            tol = self.settings['pnp_tolerance']
            res = (resP < tol and resA < tol and resB < tol)
            print(itr, resP, resA, resB)
            if not res:
                # poisson eq
                phi_old = p_alg[p_alg.settings['quantity']].copy()
                self._p_run_reactive(x=phi_old)
                phi_new = p_alg[p_alg.settings['quantity']].copy()
                # residual
                resP = np.sum(np.absolute(phi_old**2 - phi_new**2))

                # electrolyteA
                cA_old = eA_alg[eA_alg.settings['quantity']].copy()
                self._pnp_run_reactive(electrolyte='electrolyteA', x=cA_old)
                cA_new = eA_alg[eA_alg.settings['quantity']].copy()
                # residual
                resA = np.sum(np.absolute(cA_old**2 - cA_new**2))

                # electrolyteB
                cB_old = eB_alg[eB_alg.settings['quantity']].copy()
                self._pnp_run_reactive(electrolyte='electrolyteB', x=cB_old)
                cB_new = eB_alg[eB_alg.settings['quantity']].copy()
                # residual
                resB = np.sum(np.absolute(cB_old**2 - cB_new**2))
            if res:
                logger.info('Solution converged: ' + str(res))
                break

    def _pnp_run_reactive(self, electrolyte, x):
        e_alg = self.settings['electrolytes'][electrolyte]['algorithm']
        if x is None:
            x = np.zeros(shape=[e_alg.Np, ], dtype=float)
        e_alg[e_alg.settings['quantity']] = x
        rlx = e_alg.settings['relaxation_quantity']
        # Reference for residual's normalization
        e_alg._A = self._build_A(electrolyte=electrolyte, force=True)
        ref = np.sum(np.absolute(e_alg._A.diagonal())) or 1
        for itr in range(int(e_alg.settings['max_iter'])):
            e_alg[e_alg.settings['quantity']] = x
            e_alg._A = self._build_A(electrolyte=electrolyte, force=True)
            e_alg._b = e_alg._build_b(force=True)
            e_alg._apply_BCs()
            e_alg._apply_sources()
            # Compute the normalized residual
            res = np.linalg.norm(e_alg.b-e_alg.A*x)/ref
            if res >= e_alg.settings['rxn_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                x_new = e_alg._solve()
                # Relaxation
                x_new = rlx*x_new + (1-rlx)*e_alg[e_alg.settings['quantity']]
                e_alg[e_alg.settings['quantity']] = x_new
                x = x_new
            if (res < e_alg.settings['rxn_tolerance'] or
                    e_alg.settings['sources'] == []):
                x_new = x
                logger.info('Solution converged: ' + str(res))
                break
        return x_new

    def _p_run_reactive(self, x):
        p_alg = self.settings['poisson_alg']
        if x is None:
            x = np.zeros(shape=[p_alg.Np, ], dtype=float)
        p_alg[p_alg.settings['quantity']] = x
        rlx = p_alg.settings['relaxation_quantity']
        # Reference for residual's normalization
        p_alg._build_A(force=True)
        ref = np.sum(np.absolute(p_alg.A.diagonal())) or 1
        for itr in range(int(p_alg.settings['max_iter'])):
            p_alg[p_alg.settings['quantity']] = x
            p_alg._build_A(force=True)
            self._update_b_poisson()
            p_alg._apply_BCs()
            p_alg._apply_sources()
            # Compute the normalized residual
            res = np.linalg.norm(p_alg.b-p_alg.A*x)/ref
            if res >= p_alg.settings['rxn_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                x_new = p_alg._solve()
                # Relaxation
                x_new = rlx*x_new + (1-rlx)*p_alg[p_alg.settings['quantity']]
                p_alg[p_alg.settings['quantity']] = x_new
                x = x_new
            if (res < p_alg.settings['rxn_tolerance'] or
                    p_alg.settings['sources'] == []):
                x_new = x
                logger.info('Solution converged: ' + str(res))
                break
        return x_new
