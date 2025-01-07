import logging
import sys

import numpy as np
from numpy.linalg import norm
# from scipy.optimize.nonlin import TerminationCondition
from tqdm.auto import tqdm

from openpnm.algorithms import Transport
from openpnm.utils import Docorator, TypedList

__all__ = ["ReactiveTransport"]


docstr = Docorator()
logger = logging.getLogger(__name__)


def maxnorm(x):
    return np.absolute(x).max()


class TerminationCondition:
    """
    Termination condition for an iteration. It is terminated if

    - |F| < f_rtol*|F_0|, AND
    - |F| < f_tol

    AND

    - |dx| < x_rtol*|x|, AND
    - |dx| < x_tol

    Notes
    -----
    This class was deprecated from Scipy, so I dug up the old source code and
    pasted it here for now.

    """
    def __init__(self, f_tol=None, f_rtol=None, x_tol=None, x_rtol=None,
                 iter=None, norm=maxnorm):

        if f_tol is None:
            f_tol = np.finfo(np.float_).eps ** (1./3)
        if f_rtol is None:
            f_rtol = np.inf
        if x_tol is None:
            x_tol = np.inf
        if x_rtol is None:
            x_rtol = np.inf

        self.x_tol = x_tol
        self.x_rtol = x_rtol
        self.f_tol = f_tol
        self.f_rtol = f_rtol

        self.norm = norm

        self.iter = iter

        self.f0_norm = None
        self.iteration = 0

    def check(self, f, x, dx):
        self.iteration += 1
        f_norm = self.norm(f)
        x_norm = self.norm(x)
        dx_norm = self.norm(dx)

        if self.f0_norm is None:
            self.f0_norm = f_norm

        if f_norm == 0:
            return 1

        if self.iter is not None:
            # backwards compatibility with SciPy 0.6.0
            return 2 * (self.iteration > self.iter)

        # NB: condition must succeed for rtol=inf even if norm == 0
        return int((f_norm <= self.f_tol
                    and f_norm/self.f_rtol <= self.f0_norm)
                   and (dx_norm <= self.x_tol
                        and dx_norm/self.x_rtol <= x_norm))


@docstr.get_sections(base="ReactiveTransportSettings", sections=["Parameters"])
@docstr.dedent
class ReactiveTransportSettings:
    r"""

    Parameters
    ----------
    %(TransportSettings.parameters)s
    sources : list
        List of source terms that have been added
    relaxation_factor : float (default = 1.0)
        A relaxation factor to control under-relaxation for the quantity
        solving for. Factor approaching 0 leads to improved stability but
        slower simulation. Factor approaching 1 gives fast simulation but
        may be unstable.
    newton_maxiter : int
        Maximum number of iterations allowed for the nonlinear solver to
        converge.
    f_rtol : float
        Relative tolerance for the solution residual
    x_rtol : float
        Relative tolerance for the solution vector

    """
    relaxation_factor = 1.0
    newton_maxiter = 5000
    f_rtol = 1e-6
    x_rtol = 1e-6


@docstr.get_sections(base="ReactiveTransport", sections=["Parameters"])
@docstr.dedent
class ReactiveTransport(Transport):
    r"""
    A subclass for steady-state simulations with (optional) source terms.

    Parameters
    ----------
    %(Transport.parameters)s

    Notes
    -----
    This subclass performs steady simulations of transport phenomena with
    reactions when source terms are added.

    """

    def __init__(self, name="react_trans_?", **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(ReactiveTransportSettings())

    def set_source(self, pores, propname, mode="add"):
        r"""
        Applies a given source term to the specified pores

        Parameters
        ----------
        pores : array_like
            The pore indices where the source term should be applied.
        propname : str
            The property name of the source term model to be applied.
        mode : str
            Controls how the sources are applied. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'add'        (default) Adds supplied source term to already
                         existing ones.
            'remove'     Deletes given source term from the specified
                         locations.
            'clear'      Deletes given source term from all locations (ignores
                         the ``pores`` argument).
            ===========  =====================================================

        Notes
        -----
        Source terms cannot be applied in pores where boundary conditions
        have already been set. Attempting to do so will result in an error
        being raised.

        """
        # If a list of modes was sent, do them each in order
        if isinstance(mode, list):
            for item in mode:
                self.set_source(pores=pores, propname=propname, mode=item)
            return
        propname = self._parse_prop(propname, "pore")
        # Check if any BC is already set in the same locations
        locs_BC = np.zeros(self.Np, dtype=bool)
        for item in self["pore.bc"].keys():
            locs_BC = np.isfinite(self[f"pore.bc.{item}"])
        if np.any(locs_BC[pores]):
            raise Exception("BCs present in given pores, can't assign source term")
        prop = "pore.source." + propname.split(".", 1)[1]
        if mode == "add":
            if prop not in self.keys():
                self[prop] = False
            self[prop][pores] = True
        elif mode == "remove":
            self[prop][pores] = False
        elif mode == "clear":
            self[prop] = False
        else:
            raise Exception(f"Unsupported mode {mode}")

    def _apply_sources(self):
        """
        Updates ``A`` and ``b``, applying source terms to specified pores.

        Notes
        -----
        Phase and physics objects are also updated before applying source
        terms to ensure that source terms values are associated with the
        current value of 'quantity'.

        """
        try:
            phase = self.project[self.settings.phase]
            for item in self["pore.source"].keys():
                # Fetch linearized values of the source term
                Ps = self["pore.source." + item]
                S1, S2 = [phase[f"pore.{item}.{Si}"] for Si in ["S1", "S2"]]
                # Modify A and b: diag(A) += -S1, b += S2
                diag = self.A.diagonal()
                diag[Ps] += -S1[Ps]
                self.A.setdiag(diag)
                self.b[Ps] += S2[Ps]
        except KeyError:
            pass

    def _run_special(self, solver, x0, verbose=None):
        r"""
        Repeatedly updates ``A``, ``b``, and the solution guess within
        according to the applied source term then calls ``_solve`` to
        solve the resulting system of linear equations.

        Stops when the max-norm of the residual drops by at least
        ``f_rtol``:

            ``norm(R_n) < norm(R_0) * f_rtol``

        AND

            ``norm(dx) < norm(x) * x_rtol``

        where R_i is the residual at ith iteration, x is the solution at
        current iteration, and dx is the change in the solution between two
        consecutive iterations. ``f_rtol`` and ``x_rtol`` are defined in
        the algorithm's settings under: ``alg.settings['f_rtol']``, and
        ``alg.settings['x_rtol']``, respectively.

        Parameters
        ----------
        x0 : ndarray
            Initial guess of the unknown variable

        """
        w = self.settings["relaxation_factor"]
        maxiter = self.settings["newton_maxiter"]
        f_rtol = self.settings["f_rtol"]
        x_rtol = self.settings["x_rtol"]
        xold = self.x
        dx = self.x - xold
        condition = TerminationCondition(
            f_tol=np.inf, f_rtol=f_rtol, x_rtol=x_rtol, norm=norm
        )

        tqdm_settings = {
            "total": 100,
            "desc": f"{self.name} : Newton iterations",
            "disable": not verbose,
            "file": sys.stdout,
            "leave": False,
        }

        with tqdm(**tqdm_settings) as pbar:
            for i in range(maxiter):
                res = self._get_residual()
                progress = self._get_progress(res)
                pbar.update(progress - pbar.n)
                is_converged = bool(condition.check(f=res, x=xold, dx=dx))
                if is_converged:
                    pbar.update(100 - pbar.n)
                    self.soln.is_converged = is_converged
                    logger.info(f"Solution converged, residual norm: {norm(res):.4e}")
                    return
                super()._run_special(solver=solver, x0=xold, w=w)
                dx = self.x - xold
                xold = self.x
                logger.info(f"Iteration #{i:<4d} | Residual norm: {norm(res):.4e}")
                self.soln.num_iter = i + 1

        self.soln.is_converged = False
        logger.warning(f"{self.name} didn't converge after {maxiter} iterations")

    def _get_progress(self, res):
        """
        Returns an approximate value for completion percent of Newton iterations.
        """
        if not hasattr(self, "_f0_norm"):
            self._f0_norm = norm(res)
        f_rtol = self.settings.f_rtol
        norm_reduction = norm(res) / self._f0_norm / f_rtol
        progress = (1 - max(np.log10(norm_reduction), 0) / np.log10(1 / f_rtol)) * 100
        return max(0, progress)

    def _update_A_and_b(self):
        r"""
        Builds/updates A, b based on the recent solution on algorithm object.
        """
        self._update_iterative_props()
        super()._update_A_and_b()
        self._apply_sources()

    def _get_residual(self, x=None):
        r"""
        Calculates solution residual based on the given ``x`` based on the
        following formula:

            ``R = A * x - b``

        """
        if x is None:
            x = self.x
        return self.A * x - self.b

    def set_BC(self, pores=None, bctype=[], bcvalues=[], mode="add"):
        msg = "Source term already present in given pores, can't assign BCs"
        # Ensure that given pores do not have source terms already set
        try:
            for item in self["pore.source"].keys():
                if np.any(self["pore.source." + item][pores]):
                    raise Exception(msg)
        except KeyError:
            pass
        # Assign BCs if above check passes
        super().set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues, mode=mode)
