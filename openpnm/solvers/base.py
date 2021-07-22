from numpy.linalg import norm


class BaseSolver:

    def __init__(self):
        pass

    def solve(self):
        raise NotImplementedError


class DirectSolver(BaseSolver): ...


class IterativeSolver(BaseSolver):

    def __init__(self, tol=1e-8, maxiter=5000):
        self.tol = tol
        self.maxiter = maxiter

    def _get_atol(self, b):
        r"""
        Returns the absolute tolerance ``atol`` that corresponds to the
        the given tolerance ``tol``.

        Notes
        -----
        ``atol`` is calculated to satisfy the following stopping criterion:
            ``norm(A*x-b)`` <= ``atol``

        """
        return norm(self.b) * self.tol

    def _get_rtol(self, x0):
        r"""
        Returns the relative tolerance ``rtol`` that corresponds to the
        the given tolerance ``tol``.

        Notes
        -----
        ``rtol`` is defined based on the following formula:
            ``rtol = residual(@x_final) / residual(@x0)``

        """
        res0 = self._get_residual(x=x0)
        atol = self._get_atol()
        rtol = atol / res0
        return rtol

    def _get_residual(self, A, b, x):
        r"""
        Calculates the residual based on the given ``x`` using:
            ``res = norm(A*x - b)``
        """
        return norm(A*x - b)
