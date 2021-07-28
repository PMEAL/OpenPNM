    def _get_solver(self):
        r"""
        Fetch solver object based on solver settings stored in settings dict.

        Notes
        -----
        The returned object can be called via ``obj.solve(A, b, x0[optional])``

        """
        # SciPy
        if self.settings['solver_family'] == 'scipy':
            def solver(A, b, atol=None, rtol=None, max_it=None, x0=None):
                r"""
                Wrapper method for scipy sparse linear solvers.
                """
                ls = getattr(scipy.sparse.linalg, self.settings['solver_type'])
                if self.settings["solver_type"] == "spsolve":
                    x = ls(A=A, b=b)
                else:
                    tol = self.settings["solver_tol"]
                    x, _ = ls(A=A, b=b, atol=atol, tol=tol, maxiter=max_it, x0=x0)
                return x
        # PETSc
        elif self.settings['solver_family'] == 'petsc':
            def solver(A, b, atol=None, rtol=None, max_it=None, x0=None):
                r"""
                Wrapper method for PETSc sparse linear solvers.
                """
                from openpnm.utils.petsc import PETScSparseLinearSolver as SLS
                temp = {"type": self.settings["solver_type"],
                        "preconditioner": self.settings["solver_preconditioner"]}
                ls = SLS(A=A, b=b, settings=temp)
                x = ls.solve(x0=x0, atol=atol, rtol=rtol, max_it=max_it)
                return x
        # PyAMG
        elif self.settings['solver_family'] == 'pyamg':
            def solver(A, b, rtol=None, max_it=None, x0=None, **kwargs):
                r"""
                Wrapper method for PyAMG sparse linear solvers.
                """
                import pyamg
                ml = pyamg.smoothed_aggregation_solver(A)
                x = ml.solve(b=b, x0=x0, tol=rtol, maxiter=max_it, accel="bicgstab")
                return x
        # PyPardiso
        elif self.settings['solver_family'] in ['pypardiso', 'pardiso']:
            def solver(A, b, **kwargs):
                r"""
                Wrapper method for PyPardiso sparse linear solver.
                """
                import pypardiso
                x = pypardiso.spsolve(A=A, b=b)
                return x
        # CuPy
        elif self.settings['solver_family'] == 'cupy':  # pragma: no cover
            try:
                import cupy
                import cupyx.scipy.sparse.linalg
                cupyx.scipy.sparse.linalg.lschol = cupyx.linalg.sparse.lschol
            except ModuleNotFoundError:
                msg = "CuPy missing. Install via: conda install -c conda-forge cupy"
                raise Exception(msg)

            def solver(A, b, atol=None, rtol=None, max_it=None, x0=None, **kwargs):
                r"""
                Wrapper method for CuPy sparse linear solvers.
                """
                x0 = x0 if x0 is None else cupy.array(x0)
                b = cupy.array(b)
                A = cupy.sparse.csr_matrix(A)
                direct = ["spsolve", "lsqr", "lschol"]
                iterative = ["cg", "gmres"]
                solver_type = self.settings["solver_type"]
                args = {"A": A, "b": b}
                if solver_type in direct + iterative:
                    ls = getattr(cupyx.scipy.sparse.linalg, solver_type)
                    if solver_type in iterative:
                        args.update({"tol": rtol, "maxiter": max_it, "x0": x0})
                else:
                    raise Exception(f"Unsupported solver type: {solver_type}")
                out = ls(**args)
                x = out[0] if isinstance(out, tuple) else out
                return cupy.asnumpy(x)
        else:
            raise Exception(f"{self.settings['solver_family']} not available.")

        return solver