import numpy as np
from scipy.sparse.csgraph import laplacian
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class PoissonNernstPlanck(ReactiveTransport):
    r"""
    A subclass of GenericTransport to solve the Poisson Nernst-Planck equations
    """
    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'hydraulic_conductance':
                   'throat.hydraulic_conductance.solvent',
                   's_scheme': 'powerlaw',
                   'tolerance': 1e-4,
                   'max_iter': 10}
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)
        self._setup_quantities_conductances()

    def setup(self, phase=None, electrolytes=None, pressure_field=None,
              potential_field=None, **kwargs):
        r"""
        """
        if phase:
            self.settings['phase'] = phase.name
        if electrolytes:
            self.settings['electrolytes'] = electrolytes
        if pressure_field:
            self.settings['pressure_field'] = pressure_field
        if potential_field:
            self.settings['potential_field'] = potential_field
        self._setup_quantities_conductances()
        super().setup(**kwargs)

    def _setup_quantities_conductances(self):
        if self.settings['pressure_field'] is not None:
            self.settings['pressure_field'].setup(quantity='pore.pressure')
            self.settings['pressure_field'].setup(
                           conductance='throat.hydraulic_conductance.solvent')
        if self.settings['potential_field'] is not None:
            self.settings['potential_field'].setup(quantity='pore.potential')
            self.settings['potential_field'].setup(
                           conductance='throat.electrical_conductance.solvent')
        if self.settings['electrolytes'] is not None:
            for e in self.settings['electrolytes']:
                e.setup(quantity='pore.concentration.'+e.name)

    def _build_A(self, electrolyte=None, force=False):
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        f_alg = self.settings['pressure_field']
        p_alg = self.settings['potential_field']
        e_alg = electrolyte
        s_dis = self.settings['s_scheme']

        # Flow
        try:
            f_alg[f_alg.settings['quantity']]
        except KeyError:
            f_alg.run()
        P = f_alg[f_alg.settings['quantity']]

        gh = phase[self.settings['hydraulic_conductance']]
        gd = phase['throat.diffusive_conductance.'+e_alg.name]
        gd = np.tile(gd, 2)

        conns = network['throat.conns']
        Qij = -gh*np.diff(P[conns], axis=1).squeeze()
        Qij = np.append(Qij, -Qij)

        phi = p_alg[p_alg.settings['quantity']]
        S = 1/phase[p_alg.settings['conductance']]
        L = (network['throat.conduit_lengths.pore1'] +
             network['throat.conduit_lengths.throat'] +
             network['throat.conduit_lengths.pore1'])
        z = phase['pore.valence.'+e_alg.name][0]
        D = phase['pore.diffusivity.'+e_alg.name][0]
        F = 96485.3329
        R = 8.3145
        T = 298

        grad_phi = -np.diff(phi[conns], axis=1).squeeze() / L
        mig = ((z*F*D*S)/(R*T)) * grad_phi
        mig = np.append(mig, -mig)

        adv_mig = Qij-mig

        Peij_adv_mig = adv_mig/gd
        Peij_adv_mig[(Peij_adv_mig < 1e-10) & (Peij_adv_mig >= 0)] = 1e-10
        Peij_adv_mig[(Peij_adv_mig > -1e-10) & (Peij_adv_mig <= 0)] = -1e-10

        adv_mig = Peij_adv_mig*gd

        if force:
            self._pure_A = None
        if self._pure_A is None:
            if s_dis == 'upwind':
                w = gd + np.maximum(0, -adv_mig)
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'hybrid':
                w = np.maximum(0, np.maximum(-adv_mig, gd-adv_mig/2))
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'powerlaw':
                w = gd * np.maximum(0, (1 - 0.1*np.abs(Peij_adv_mig))**5) + \
                    np.maximum(0, -adv_mig)
                A = network.create_adjacency_matrix(weights=w)
            elif s_dis == 'exponential':
                w = -adv_mig / (1 - np.exp(Peij_adv_mig))
                A = network.create_adjacency_matrix(weights=w)
            else:
                raise Exception('Unrecognized discretization scheme: ' + s_dis)

            A = laplacian(A)
            self._pure_A = A
        self.A = self._pure_A.copy()
        return self.A

    def _update_b_poisson(self):
        phase = self.project.phases()[self.settings['phase']]
        p_alg = self.settings['potential_field']
        e_alg = self.settings['electrolytes']

        Sum = 0
        for e in e_alg:
            Sum += phase['pore.valence.'+e.name] * e[e.settings['quantity']]
        p_alg._b = Sum.astype('float64')

    def run(self):
        p_alg = self.settings['potential_field']
        e_alg = self.settings['electrolytes']

        p_alg[p_alg.settings['quantity']] = 0
        for e in e_alg:
            e[e.settings['quantity']] = 0

        tol = self.settings['tolerance']
        # Initialize the residuals
        res = {}
        res['potential'] = 1e+06
        for e in e_alg:
            res[e.name] = 1e+06

        for itr in range(int(self.settings['max_iter'])):
            print('\n', itr, end='\t')
            r = str([float(format(i, '.3g')) for i in res.values()])[1:-1]
            logger.info('Iter: ' + str(itr) + ', Tols: ' + r)
            convergence = max(i for i in res.values()) < tol
            if not convergence:
                # poisson eq
                phi_old = p_alg[p_alg.settings['quantity']].copy()
                self._run_reactive(alg=p_alg, x=phi_old)
                phi_new = p_alg[p_alg.settings['quantity']].copy()
                # residual
                res['potential'] = np.sum(np.absolute(phi_old**2 - phi_new**2))

                # electrolytes
                for e in e_alg:
                    c_old = e[e.settings['quantity']].copy()
                    self._run_reactive(alg=e, x=c_old)
                    c_new = e[e.settings['quantity']].copy()
                    # residual
                    res[e.name] = np.sum(np.absolute(c_old**2 - c_new**2))

            if convergence:
                logger.info('Solution converged: ' + str(res))
                break

    def _run_reactive(self, alg, x):
        if x is None:
            x = np.zeros(shape=[alg.Np, ], dtype=float)
        alg[alg.settings['quantity']] = x
        rlx = alg.settings['relaxation_quantity']

        # Reference for residual's normalization
        if alg == self.settings['potential_field']:
            alg._build_A(force=True)
        else:
            alg._A = self._build_A(electrolyte=alg, force=True)
        ref = np.sum(np.absolute(alg._A.diagonal())) or 1

        for itr in range(int(alg.settings['max_iter'])):
            alg[alg.settings['quantity']] = x

            if alg == self.settings['potential_field']:
                alg._build_A(force=True)
                self._update_b_poisson()
            else:
                alg._A = self._build_A(electrolyte=alg, force=True)
                alg._b = alg._build_b(force=True)

            alg._apply_BCs()
            alg._apply_sources()
            # Compute the normalized residual
            res = np.linalg.norm(alg.b-alg.A*x)/ref
            if res >= alg.settings['rxn_tolerance']:
                logger.info('Tolerance not met: ' + str(res))
                x_new = alg._solve()
                # Relaxation
                x_new = rlx*x_new + (1-rlx)*alg[alg.settings['quantity']]
                alg[alg.settings['quantity']] = x_new
                x = x_new
            if (res < alg.settings['rxn_tolerance'] or
                    alg.settings['sources'] == []):
                x_new = x
                logger.info('Solution converged: ' + str(res))
                break
        return x_new
