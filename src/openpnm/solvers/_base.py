from numpy.linalg import norm

__all__ = ['BaseSolver', 'DirectSolver', 'IterativeSolver']


class BaseSolver:
    """Base class for all solvers."""
    def __init__(self):
        ...

    def solve(self):
        """Solves the given linear system of equations Ax=b."""
        raise NotImplementedError


class DirectSolver(BaseSolver):
    """Base class for all direct solvers."""
    ...


class IterativeSolver(BaseSolver):
    """Base class for iterative solvers."""
    def __init__(self, tol=1e-8, maxiter=1000):
        self.tol = tol
        self.maxiter = maxiter
        self.atol = None        # needs to be evaluated later
        self.rtol = None        # needs to be evaluated later

    def _get_atol(self, b):
        r"""
        Returns the absolute tolerance ``atol`` that corresponds to the
        the given tolerance ``tol``.

        Notes
        -----
        ``atol`` is calculated to satisfy the following stopping criterion:
            ``norm(A*x-b)`` <= ``atol``

        """
        return norm(b) * self.tol

    def _get_rtol(self, A, b, x0):
        r"""
        Returns the relative tolerance ``rtol`` that corresponds to the
        the given tolerance ``tol``.

        Notes
        -----
        ``rtol`` is defined based on the following formula:
            ``rtol = residual(@x_final) / residual(@x0)``

        """
        res0 = self._get_residual(A, b, x0)
        atol = self._get_atol(b)
        rtol = atol / res0
        return rtol

    def _get_residual(self, A, b, x):
        r"""
        Calculates the residual based on the given ``x`` using:
            ``res = norm(A*x - b)``
        """
        return norm(A * x - b)
