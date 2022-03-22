import logging
import sys
import numpy as np
import scipy as sp
import scipy.sparse
from openpnm.solvers import IterativeSolver
logger = logging.getLogger(__name__)
try:
    import petsc4py
    # Next line must be before importing PETSc
    petsc4py.init(sys.argv)
    from petsc4py import PETSc
except ModuleNotFoundError:
    pass


__all__ = ['PETScLinearSolver']


class PETScLinearSolver(IterativeSolver):
    r"""
    Solves the sparse linear system Ax = b using petsc solvers.

    Notes
    -----
    Parallel computing is supported and matrix partitioning over the
    available cores is automatically handled by running:

    .. code::

        $ mpirun -np num_cores python script.py

    where ``num_cores`` must be substituted with the number of cores.

    """
    def _create_solver(self):
        r"""
        This method creates the petsc sparse linear solver.
        """
        # https://petsc.org/release/docs/manualpages/KSP/KSPType.html
        iterative = [
            'richardson', 'chebyshev', 'cg', 'groppcg', 'pipecg', 'pipecgrr',
            'cgne', 'nash', 'stcg', 'gltr', 'fcg', 'pipefcg', 'gmres',
            'pipefgmres', 'fgmres', 'lgmres', 'dgmres', 'pgmres', 'tcqmr',
            'bcgs', 'ibcgs', 'fbcgs', 'fbcgsr', 'bcgsl', 'pipebcgs', 'cgs',
            'tfqmr', 'cr', 'pipecr', 'lsqr', 'preonly', 'qcg', 'bicg',
            'minres', 'symmlq', 'lcd', 'python', 'gcr', 'pipegcr', 'tsirm',
            'cgls', 'fetidp']
        # https://petsc.org/release/docs/manualpages/PC/PCType.html
        preconditioners = [
            'none', 'jacobi', 'sor', 'lu', 'shell', 'bjacobi', 'mg',
            'eisenstat', 'ilu', 'icc', 'asm', 'gasm', 'ksp', 'composite',
            'redundant', 'spai', 'nn', 'cholesky', 'pbjacobi', 'mat', 'hypre',
            'parms', 'fieldsplit', 'tfs', 'ml', 'galerkin', 'exotic', 'cp',
            'bfbt', 'lsc', 'python', 'pfmg', 'syspfmg', 'redistribute', 'svd',
            'gamg', 'sacusp', 'sacusppoly', 'bicgstabcusp', 'ainvcusp',
            'chowiluviennacl', 'rowscalingviennacl', 'saviennacl', 'bddc',
            'kaczmarz', 'telescope']
        direct_lu = ['mumps', 'superlu_dist', 'umfpack', 'klu']
        direct_cholesky = ['mumps', 'cholmod']
        valid_solvers = iterative + direct_lu + direct_cholesky

        solver = self.solver_type
        preconditioner = self.preconditioner

        if solver not in valid_solvers:
            raise Exception(f"{solver} solver not availabe, choose another solver")
        if preconditioner not in preconditioners:
            raise Exception(f"{preconditioner} not found, choose another preconditioner")

        self.ksp = PETSc.KSP()
        self.ksp.create(PETSc.COMM_WORLD)

        if solver in direct_lu:
            self.ksp.getPC().setType('lu')
            self.ksp.getPC().setFactorSolverType(solver)
            self.ksp.setType('preonly')
        elif solver in direct_cholesky:
            self.ksp.getPC().setType('cholesky')
            self.ksp.getPC().setFactorSolverType(solver)
            self.ksp.setType('preonly')
        elif solver in preconditioners:
            self.ksp.getPC().setType(solver)
            self.ksp.setType('preonly')
        elif solver in iterative:
            self.ksp.getPC().setType(preconditioner)
            self.ksp.setType(solver)

    def _set_tolerances(self, atol=None, rtol=None, maxiter=None):
        r"""
        Set absolute and relative tolerances, and maximum number of iterations.
        """
        atol = self.atol if atol is None else atol
        rtol = self.rtol if rtol is None else rtol
        maxiter = self.maxiter if maxiter is None else maxiter
        # BUG: PETSc misses rtol requirement by ~10-20X -> Report to petsc4py
        self.ksp.setTolerances(atol=None, rtol=rtol/50, max_it=maxiter)

    def _assemble_A(self):
        r"""
        This method creates the petsc sparse coefficients matrix from the
        OpenPNM scipy one. The method also equally decomposes the matrix at
        certain rows into different blocks (each block contains all the
        columns) and distributes them over the pre-assigned cores for parallel
        computing. The method can be used in serial.
        """

        # Create a petsc sparse matrix
        self.petsc_A = PETSc.Mat()
        self.petsc_A.create(PETSc.COMM_WORLD)
        self.petsc_A.setSizes([self.m, self.n])
        self.petsc_A.setType('aij')  # sparse
        self.petsc_A.setUp()

        # Loop over owned block of rows on this processor
        # and insert entry values (for parallel computing).
        self.Istart, self.Iend = self.petsc_A.getOwnershipRange()

        # Assign values to the coefficients matrix from the scipy sparse csr
        size_tmp = self.A.shape
        # Row indices
        csr1 = self.A.indptr[self.Istart:self.Iend+1] - self.A.indptr[self.Istart]
        ind1 = self.A.indptr[self.Istart]
        ind2 = self.A.indptr[self.Iend]
        csr2 = self.A.indices[ind1:ind2]  # column indices
        csr3 = self.A.data[ind1:ind2]  # data
        self.petsc_A = PETSc.Mat().createAIJ(size=size_tmp,
                                             csr=(csr1, csr2, csr3))
        # Communicate off-processor values and setup internal data structures
        # for performing parallel operations
        self.petsc_A.assemblyBegin()
        self.petsc_A.assemblyEnd()

    def _assemble_b_and_x(self):
        r"""
        Initialize the solution vector (self.petsc_x), which is a dense
        matrix (1D vector) and defines the rhs vector (self.petsc_b) from
        the existing data.
        """
        # Get vector(s) compatible with the matrix (same parallel layout)
        # passing same communicator as the A matrix

        # Global solution vector (all the local solutions will return to it)
        self.petsc_s = PETSc.Vec()
        self.petsc_s.create(PETSc.COMM_WORLD)
        self.petsc_s.setSizes(self.m)
        self.petsc_s.setFromOptions()

        self.Istart, self.Iend = self.petsc_s.getOwnershipRange()

        self.petsc_x = (self.petsc_s).duplicate()
        self.petsc_b = (self.petsc_s).duplicate()

        # Set the solution vector to zeros or the given initial guess (if any)
        PETSc.Vec.setArray(self.petsc_x, self.x0[self.Istart: self.Iend])

        # Define the petsc rhs vector from the numpy one
        PETSc.Vec.setArray(self.petsc_b, self.b[self.Istart: self.Iend])

    def solve(self, A, b, x0=None, solver_type='cg', precondioner='jacobi',
              maxiter=None, atol=None, rtol=None):
        r"""
        Solves and returns the solution to the linear system, Ax = b.

        This method converts the solution vector from a PETSc.Vec
        instance to a numpy array, and finally destroys all the PETSc
        objects to free memory.

        Parameters
        ----------
        A : csr_matrix
            Coefficients matrix in Ax = b
        b : ndarray
            Right-hand-side vector in Ax = b
        solver_type : str, optional
            Default is the iterative solver 'cg' based on the
            Conjugate Gradient method.
        preconditioner: str, optional
            Default is the 'jacobi' preconditioner, i.e., diagonal
            scaling preconditioning. The preconditioner is used with
            iterative solvers. When a direct solver is used, this
            parameter is ignored.
        factorization_type : str, optional
            The factorization type used with the direct solver.
            Default is 'lu'. This parameter is ignored when an
            iterative solver is used.

        Returns
        -------
        ndarray
            The solution to Ax = b

        Notes
        -----
        Certain combinations of iterative solvers and precondioners
        or direct solvers and factorization types are not supported.
        The summary table of the different possibilities
        can be found
        `here <https://petsc.org/main/overview/linear_solve_table>`_

        """
        self.b = b
        self.A = sp.sparse.csr_matrix(A)
        self.m, self.n = self.A.shape
        self.x0 = np.zeros_like(self.b) if x0 is None else x0

        self.solver_type = solver_type
        self.preconditioner = precondioner

        self.atol = self._get_atol(self.b)
        self.rtol = self._get_rtol(self.x0)

        self._assemble_b_and_x()
        self._assemble_A()

        self._create_solver()
        self._set_tolerances(atol=atol, rtol=rtol, maxiter=maxiter)
        self.ksp.setOperators(self.petsc_A)
        self.ksp.setFromOptions()

        # Solve the linear system
        self.ksp.solve(self.petsc_b, self.petsc_x)

        # Gather the solution to all processors
        gather_to_0, self.petsc_s = PETSc.Scatter().toAll(self.petsc_x)
        gather_to_0.scatter(self.petsc_x, self.petsc_s,
                            PETSc.InsertMode.INSERT, PETSc.ScatterMode.FORWARD)

        # Convert solution vector from PETSc.Vec instance to a numpy array
        self.solution = PETSc.Vec.getArray(self.petsc_s)

        # Destroy petsc solver, coefficients matrix, rhs, and solution vectors
        PETSc.KSP.destroy(self.ksp)
        PETSc.Mat.destroy(self.petsc_A)
        PETSc.Vec.destroy(self.petsc_b)
        PETSc.Vec.destroy(self.petsc_x)
        PETSc.Vec.destroy(self.petsc_s)

        # FIXME: fetch exit_code somehow from petsc
        exit_code = 0

        return self.solution, exit_code
