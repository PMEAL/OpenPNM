import sys
import logging
import numpy as np
from numpy.linalg import norm
from scipy.optimize.nonlin import TerminationCondition
from openpnm.algorithms import GenericTransport
from openpnm.utils import TypedList, Docorator
from tqdm import tqdm


__all__ = ['ReactiveTransport']


docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='ReactiveTransportSettings', sections=['Parameters'])
@docstr.dedent
class ReactiveTransportSettings:
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
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
    sources = TypedList(types=[str])
    newton_maxiter = 5000
    f_rtol = 1e-6
    x_rtol = 1e-6


@docstr.get_sections(base='ReactiveTransport', sections=['Parameters'])
@docstr.dedent
class ReactiveTransport(GenericTransport):
    r"""
    A subclass for steady-state simulations with (optional) source terms.

    Parameters
    ----------
    %(GenericTransport.parameters)s

    Notes
    -----
    This subclass performs steady simulations of transport phenomena with
    reactions when source terms are added.

    """

    def __init__(self, name='react_trans_#', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(ReactiveTransportSettings())

    def set_source(self, propname, pores, mode='overwrite'):
        r"""
        Applies a given source term to the specified pores

        Parameters
        ----------
        propname : str
            The property name of the source term model to be applied.
        pores : array_like
            The pore indices where the source term should be applied.
        mode : str
            Controls how the sources are applied (see table under Notes).
            The default is 'overwrite'. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'add'        Adds supplied source term to already existing
                         ones.
            'overwrite'  Deletes all existing source terms of the given
                         ``propname`` then adds the specified new ones.
            ===========  =====================================================

        Notes
        -----
        Source terms cannot be applied in pores where boundary conditions
        have already been set. Attempting to do so will result in an error
        being raised.

        """
        propname = self._parse_prop(propname, "pore")
        locs = self.to_mask(pores=pores)
        # Check if any BC is already set in the same locations
        locs_BC = np.isfinite(self['pore.bc_value']) + np.isfinite(self['pore.bc_rate'])
        if (locs & locs_BC).any():
            raise Exception("BCs present in given pores, can't assign source term")
        if mode == 'overwrite':
            self[propname] = False
        if mode == 'add':
            if propname not in self.keys():
                self[propname] = False
        self[propname][locs] = True
        # Check if propname already in source term list
        if propname not in self.settings['sources']:
            self.settings['sources'].append(propname)

    def _apply_sources(self):
        """
        Updates ``A`` and ``b``, applying source terms to specified pores.

        Notes
        -----
        Phase and physics objects are also updated before applying source
        terms to ensure that source terms values are associated with the
        current value of 'quantity'.

        """
        phase = self.project[self.settings.phase]
        for item in self.settings['sources']:
            # Fetch linearized values of the source term
            Ps = self[item]
            S1, S2 = [phase[f"{item}.{Si}"] for Si in ["S1", "S2"]]
            # Modify A and b: diag(A) += -S1, b += S2
            diag = self.A.diagonal()
            diag[Ps] += -S1[Ps]
            self.A.setdiag(diag)
            self.b[Ps] += S2[Ps]

    def _run_special(self, solver, x0, verbose=True):
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
        w = self.settings['relaxation_factor']
        maxiter = self.settings['newton_maxiter']
        f_rtol = self.settings['f_rtol']
        x_rtol = self.settings['x_rtol']
        xold = self.x
        dx = self.x - xold
        condition = TerminationCondition(f_rtol=f_rtol, x_rtol=x_rtol)

        tqdm_settings = {
            "total": 100,
            "desc": f"{self.name} : Newton iterations",
            "disable": not verbose,
            "file": sys.stdout,
            "leave": False
        }

        with tqdm(**tqdm_settings) as pbar:
            for i in range(maxiter):
                self.soln.num_iter = i + 1
                res = self._get_residual()
                progress = self._get_progress(res)
                pbar.update(progress - pbar.n)
                is_converged = bool(condition.check(f=res, x=xold, dx=dx))
                if is_converged:
                    pbar.update(100 - pbar.n)
                    self.soln.is_converged = is_converged
                    logger.info(f'Solution converged, residual norm: {norm(res):.4e}')
                    return
                super()._run_special(solver=solver, x0=xold, w=w)
                dx = self.x - xold
                xold = self.x
                logger.info(f'Iteration #{i:<4d} | Residual norm: {norm(res):.4e}')

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
        progress = (1 - max(np.log10(norm_reduction), 0) / np.log10(1/f_rtol)) * 100
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

    @docstr.dedent
    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        Applies boundary conditions to specified pores if no source terms
        are already assigned to these pores. Otherwise, raise an error.

        Parameters
        ----------
        %(GenericTransport._set_BC.parameters)s

        Notes
        -----
        %(GenericTransport._set_BC.notes)s

        """
        msg = "Source term already present in given pores, can't assign BCs"
        # Ensure that given pores do not have source terms already set
        for item in self.settings['sources']:
            if np.any(self[item][pores]):
                raise Exception(msg)
        # Assign BCs if above check passes
        super()._set_BC(pores=pores, bctype=bctype, bcvalues=bcvalues, mode=mode)
