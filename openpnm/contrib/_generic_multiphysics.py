import logging
import numpy as np
from openpnm.utils import Docorator
from openpnm.algorithms import Algorithm
from openpnm.contrib import TransientMultiPhysics
from openpnm.algorithms._solution import SolutionContainer, TransientSolution


logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = ['GenericMultiPhysics']


@docstr.dedent
class GenericMultiPhysicsSettings:
    r"""

    Parameters
    ----------
    transient_algorithms: list
        List of transient algorithm objects to be solved in a coupled manner
    steady_algorithms: list
        List of steady state algorithm objects that depend on transient
        solution

    """
    transient_algorithms = []
    steady_algorithms = []


@docstr.dedent
class GenericMultiPhysics(Algorithm):
    r"""
    A subclass for generic multiphysics simulations.  By generic, we mean the
    algorithms passed can be either transient or steady state algorithms. In
    some cases, the steady state algorithm can depend on the solution to the
    transient algorithm and therefore, needs to be re-evaluated at each time
    iteration. This generic multiphysics class enables simulations of this type
    coupling trasient and non-transient algorithms whose physics depend on one
    another.

    """

    def __init__(self, transient_algorithms, steady_algorithms, **kwargs):
        super().__init__(**kwargs)
        self.settings._update(GenericMultiPhysicsSettings())
        self.settings.transient_algorithms = [alg.name for alg in 
                                              transient_algorithms]
        self.settings.steady_algorithms = [alg.name for alg in 
                                           steady_algorithms]
        self._transient_algs = transient_algorithms
        self._steady_algs = steady_algorithms
 

    def run(self, x0, tspan, interval=None, integrator=None, g_tol=None,
            g_max_iter=None):
        """
        Runs all of the transient and steady state algorithms in a coupled
        manner. Steady state algorithms are updated after each time iteration
        which is one gummel iteration.

        Parameters
        ----------
        x0 : ndarray or float
            Array (or scalar) containing initial condition values. For steady
            algorithms, the initial condition gets used as the initial guess.
            It is assumed that the initial condition for transient algorithms
            is first followed by the initial guess for steady state algorithms.
        tspan : array_like
            Tuple (or array) containing the integration time span.
        interval : float, optional
            The interval at which the solution is to be stored.
        integrator : Integrator, optional
            Integrator object which will be used to do the time stepping.
            Can be instantiated using openpnm.integrators module.
        g_tol : float (default = 1e-4)
            The tolerance to use for stopping Gummel iterations
        g_max_iter : int (default = 10)
            The maximum number of times to perform the Gummel iteration

        Returns
        -------

        """
        logger.info('Running GenericMultiphysics')
        network = self.project.network
        transient_algs = self._transient_algs
        steady_algs = self._steady_algs
        algs = transient_algs + steady_algs
        # split up x0
        nt = len(transient_algs)
        transient_x0 = x0[:nt*network.Np]
        steady_x0 = x0[nt*network.Np:]
        # define transient multiphysics object for transient algs
        tmp = TransientMultiPhysics(algorithms=transient_algs, network=network)
        # define array of times to iterate through
        times = np.arange(tspan[0]+interval, tspan[1], interval)
        times = np.append(times, tspan[1])  # ensure end time is included
        # initialize starting time
        t0 = tspan[0]
        # initialize solution container
        self.soln = SolutionContainer()
        for alg in algs:
            # transient solution for each alg
            x0 = alg.x.reshape((alg.Np, 1))
            print(t0)
            soln = TransientSolution(t0, x0)
            self.soln[alg.settings['quantity']] = soln
        phase = network.project.phases[0]
        for t in times:
            print(t0)
            # Initialize residuals & old/new fields for time marching
            g_res = {}
            g_old = {}
            g_new = {}
            for alg in algs:
                g_res[alg.name] = 1e6
                g_old[alg.name] = None
                g_new[alg.name] = None
            # for gummel iteration
            for itr in range(g_max_iter):
                g_r = [float(format(i, '.3g')) for i in g_res.values()]
                g_r = str(g_r)[1:-1]
                print(f'Start Gummel iter: {str(itr+1)} residuals: {g_r}')
                g_convergence = max(i for i in g_res.values()) < g_tol
                if not g_convergence:
                    # save current solution as old
                    for alg in algs:  # Save the current fields
                        g_old[alg.name] = alg.x.copy()
                    # run transient algorithms simultaneously
                    tmp.run(transient_x0, tspan=(t0, t), saveat=t,
                            integrator=integrator)
                    # run each steady algorithm consecutively
                    for i, alg in enumerate(steady_algs):
                        alg.run(x0=steady_x0[i*alg.Np:(i+1)*alg.Np])
                    # save new solution and calculate residual
                    for alg in algs:  # Save the current fields
                        g_new[alg.name] = alg.x.copy()
                        # calculate residual
                        res = g_old[alg.name] - g_new[alg.name]
                        g_res[alg.name] = np.sum(res**2)
                # check for convergence
                if g_convergence:
                    print(f'Solution for time step: {str(t)} s converged')
                    break
                if itr + 1 == g_max_iter:
                    print('Maximum number of Gummel iterations reached')
            # update initial condition after each time step
            for i, alg in enumerate(transient_algs):
                transient_x0[i*alg.Np:(i+1)*alg.Np] = alg.x
            for i, alg in enumerate(steady_algs):
                steady_x0[i*alg.Np:(i+1)*alg.Np] = alg.x
            # update t0
            t0 = t
            # save solution after each time step
            for alg in algs:
                # stack old and new solutions
                x_old = self.soln[alg.settings['quantity']]
                x_new = alg.x.reshape((alg.Np, 1))
                xs = np.hstack((x_old, x_new))
                # append new time to previous time steps
                ts = np.append(x_old.t, t)
                # write new TransientSolution
                soln = TransientSolution(ts, xs)
                self.soln[alg.settings['quantity']] = soln
