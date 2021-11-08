import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sectionsf('MultiphysicsSolverSettings',
                      sections=['Parameters'])
@docstr.dedent
class MultiphysicsSolverSettings(GenericSettings):
    r"""
    The following decribes the settings associated with the
    MultiphysicsSolver.

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
    algorithms = []
    g_tol = 1e-8
    g_max_iter = 10


class MultiphysicsSolver(GenericAlgorithm):
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
        c = MultiphysicsSolverSettings()
        self.settings._update_settings_and_docs(c)
        settings['phase'] = phase.name
        self.settings.update(settings)

    @docstr.dedent
    def setup(self, phase=None, algorithms=[], g_tol=None, g_max_iter=None, 
              **kwargs):
        r"""

        Parameters
        ----------
        %(MultiphysicsSolverSettings.parameters)s
        """
        if phase:
            self.settings['phase'] = phase.name
        if algorithms:
            self.settings['algorithms'] = algorithms
        if g_tol:
            self.settings['g_tol'] = g_tol
        if g_max_iter:
            self.settings['g_max_iter'] = g_max_iter

    def _get_algorithms(self):
        r"""
        
        Helper method to find algorithm object from algorithm names stored in
        settings. Returns list of algorithm objects.
        """
        algorithms = self.settings['algorithms']
        algs = [self.project.algorithms()[algorithms[i]] for i in
                range(len(algorithms))]
        return algs

    def _set_cache_to_false(self, algorithms):
        r"""
        
        Helper method used to set cache_A and cache_B to false in algorithm
        settings.
        """
        algs = algorithms
        for alg in algs:
            alg.settings.update({'cache_A': False, 'cache_b': False}) 
            
    def _set_quantity_to_algorithm(self, phase, algorithms):
        r"""
        
        Helper method used to set quantity to algorithm object if not already
        defined by the user. Used for setting initial conditions to algorithm
        object. If initial conditions are not specified on phase object then 
        zero is assumed.
        """
        algs = algorithms
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
            
    def run(self, t=None):
        r"""
        """
        print('―'*80)
        print('Running IonicTransport')
        # Phase, physics, and algorithms
        phase = self.project.phases()[self.settings['phase']]
        phys = self.project.find_physics(phase=phase)
        algs = self._get_algorithms()      
        # Force cache_A and cache_B to False
        self._set_cache_to_false(algorithms=algs)        
        # Define initial conditions (if not defined by the user)
        self._set_quantity_to_algorithm(phase=phase, algorithms=algs)
        # Initialize residuals & old/new fields for Gummel iterations
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
                # Solve each algorithm consecutively
                for alg in algs:
                    g_old[alg.name] = (alg[alg.settings['quantity']].copy())
                    alg._run_reactive(x0=g_old[alg.name])
                    g_new[alg.name] = (alg[alg.settings['quantity']].copy())
                    # Residual
                    g_res[alg.name] = np.sum(np.absolute(
                        g_old[alg.name]**2-g_new[alg.name]**2))
                    # update phase and regenerate physics
                    phase.update(alg.results())
                    for obj in phys:
                        obj.regenerate_models()

            if g_convergence:
                print('Solution converged')
                break
