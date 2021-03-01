import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='NernstPlanckMultiphysicsSolverSettings',
                     sections=['Parameters'])
@docstr.dedent
class NernstPlanckMultiphysicsSolverSettings(GenericSettings):
    r"""
    The following decribes the settings associated with the
    NernstPlanckMultiphysicsSolver.

    Parameters
    ----------
    phase : OpenPNM Phase object
        The object representing the solvent phase
    potential_field : str
        The name of the potential or voltage field data that is being solved
        for
    ions : list of OpenPNM object names
        The name of the ionic species being solved for
    g_tol : float (default = 1e-4)
        The tolerance to use for stopping Gummel iterations
    g_max_iter : int (default = 10)
        The maximum number if times to perform the Gummel iteration

    """
    phase = None
    potential_field = ''
    ions = []
    g_tol = 1e-8
    g_max_iter = 10


class NernstPlanckMultiphysicsSolver(GenericAlgorithm):
    r"""
    A multiphysics solver to solve the Nernst-Planck and Ionic Conduction
    system.

    Warnings
    --------
    This is not a true OpenPNM algorithm. This solver wraps the provided
    Nernst-Planck and ionic conduction algorithms and solves the associated
    system of equations.
    """

    def __init__(self, phase=None, settings={},  **kwargs):
        super().__init__(**kwargs)
        c = NernstPlanckMultiphysicsSolverSettings()
        self.settings._update_settings_and_docs(c)
        settings['phase'] = phase.name
        self.settings.update(settings)

    @docstr.dedent
    def setup(self, phase=None, potential_field='', ions=[], g_tol=None,
              g_max_iter=None, **kwargs):
        r"""

        Parameters
        ----------
        %(NernstPlanckMultiphysicsSolverSettings.parameters)s
        """
        if phase:
            self.settings['phase'] = phase.name
        if potential_field:
            self.settings['potential_field'] = potential_field
        if ions:
            self.settings['ions'] = ions
        if g_tol:
            self.settings['g_tol'] = g_tol
        if g_max_iter:
            self.settings['g_max_iter'] = g_max_iter

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
            alg.settings.update({'cache_A': False, 'cache_b': False})
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
        g_tol = self.settings['g_tol']
        g_res = {}
        g_old = {}
        g_new = {}
        for alg in algs:
            g_res[alg.name] = 1e+06
            g_old[alg.name] = None
            g_new[alg.name] = None

        # Iterate (Gummel) until solutions converge
        for itr in range(int(self.settings['g_max_iter'])):
            g_r = [float(format(i, '.3g')) for i in g_res.values()]
            g_r = str(g_r)[1:-1]
            print('Gummel iter: '+str(itr+1)+', residuals: '+g_r)
            g_convergence = max(i for i in g_res.values()) < g_tol
            if not g_convergence:
                # Ions
                for e in e_alg:
                    g_old[e.name] = (e[e.settings['quantity']].copy())
                    e._run_reactive(x0=g_old[e.name])
                    g_new[e.name] = (e[e.settings['quantity']].copy())
                    # Residual
                    g_res[e.name] = np.sum(np.absolute(
                        g_old[e.name]**2-g_new[e.name]**2))
                    phase.update(e.results())

                # Poisson eq
                for obj in phys:
                    obj.regenerate_models()
                g_old[p_alg.name] = p_alg[p_alg.settings['quantity']].copy()
                p_alg._run_reactive(x0=g_old[p_alg.name])
                g_new[p_alg.name] = p_alg[p_alg.settings['quantity']].copy()
                # Residual
                g_res[p_alg.name] = np.sum(np.absolute(
                    g_old[p_alg.name]**2 - g_new[p_alg.name]**2))
                # Update phase and physics
                phase.update(p_alg.results())
                for obj in phys:
                    obj.regenerate_models()

            if g_convergence:
                print('Solution converged')
                break
