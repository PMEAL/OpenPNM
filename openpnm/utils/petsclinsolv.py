# -*- coding: utf-8 -*-
"""
===============================================================================
openpnm_petsc_sparse_linear_solvers: A class for solving sparse linear systems
===============================================================================

"""
import scipy as sp
import sys
# petsc imports
import petsc4py
from petsc4py import PETSc
petsc4py.init(sys.argv)


class petscSparseLinearSolver:
    r"""
    Solve the sparse linear system Ax = b using petsc solvers. Parallel
    computing is supported and matrix partitioning over the available cores is
    automatically handled by running:
    $ mpirun -np 4 python3.5 script.py
    for parallel computing.
    """
    def __init__(self, A, b):
        r"""
        Initialize the sparse system of linear equations.

        Parameters
        ----------
        CoefMat : sparse matrix
            2D Coefficient matrix
        rhs : dense matrix
            1D RHS vector
        """
        self.A = A
        self.b = b

    def _initialize_petsc_sparse_coef_mat(self):
        r"""
        This method creates the petsc sparse coefficients matrix from the
        OpenPNM scipy one. The method also equally decomposes the matrix at
        certain rows into different blocks (each block contains all the
        columns) and distributes them over the pre-assigned cores for parallel
        computing. The method can be used in serial.
        """
        # Matrix of coefficients size
        self.m, self.n = (self.A).shape

        # Create a petsc sparse matrix
        self.petsc_A = PETSc.Mat()
        self.petsc_A.create(PETSc.COMM_WORLD)
        self.petsc_A.setSizes([self.m, self.n])
        self.petsc_A.setType('aij')  # sparse
        self.petsc_A.setUp()

        # Pre-allocate memory for the coefficients matrix.
        # Optional, but needed in the case where the
        # matrix does not already exist.
        # self.petsc_A.setPreallocationNNZ([ \
        # sp.sparse.csr_matrix.getnnz(A,axis=0)])

        # Loop over owned block of rows on this processor
        # and insert entry values (for parallel computing).
        self.Istart, self.Iend = self.petsc_A.getOwnershipRange()

        # Assign values to the coefficients matrix from the scipy
        # sparse csr one: petscMat = \
        # PETSc.Mat().createAIJ(size=existingMat.shape, csr= \
        # (existingMat.indptr,existingMat.indices,existingMat.data))

        size_tmp = self.A.shape
        csr1 = (self.A.indptr[self.Istart:self.Iend+1] -
                self.A.indptr[self.Istart])
        ind1 = self.A.indptr[self.Istart]
        ind2 = self.A.indptr[self.Iend]
        csr2 = self.A.indices[ind1:ind2]
        csr3 = self.A.data[ind1:ind2]
        self.petsc_A = PETSc.Mat().createAIJ(size=size_tmp,
                                             csr=(csr1, csr2, csr3))

        # Communicate off-processor values
        # and setup internal data structures
        # for performing parallel operations

        self.petsc_A.assemblyBegin()
        self.petsc_A.assemblyEnd()

    def _create_petsc_sparse_linear_solver(self):
        r"""
        This method creates the petsc sparse linear solver.
        """
        # http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
        petsc_iterative_solvers = ['richardson', 'chebyshev', 'cg', 'groppcg',
                                   'pipecg', 'pipecgrr', 'cgne', 'nash',
                                   'stcg', 'gltr', 'fcg', 'pipefcg', 'gmres',
                                   'pipefgmres', 'fgmres', 'lgmres', 'dgmres',
                                   'pgmres', 'tcqmr', 'bcgs', 'ibcgs', 'fbcgs',
                                   'fbcgsr', 'bcgsl', 'pipebcgs', 'cgs',
                                   'tfqmr', 'cr', 'pipecr', 'lsqr', 'preonly',
                                   'qcg', 'bicg', 'minres', 'symmlq', 'lcd',
                                   'python', 'gcr', 'pipegcr', 'tsirm', 'cgls',
                                   'fetidp']

        # http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType
        petsc_preconditioner_methods = ['none', 'jacobi', 'sor', 'lu', 'shell',
                                        'bjacobi', 'mg', 'eisenstat', 'ilu',
                                        'icc', 'asm', 'gasm', 'ksp',
                                        'composite', 'redundant', 'spai', 'nn',
                                        'cholesky', 'pbjacobi', 'mat', 'hypre',
                                        'parms', 'fieldsplit', 'tfs', 'ml',
                                        'galerkin', 'exotic', 'cp', 'bfbt',
                                        'lsc', 'python', 'pfmg', 'syspfmg',
                                        'redistribute', 'svd', 'gamg',
                                        'sacusp', 'sacusppoly', 'bicgstabcusp',
                                        'ainvcusp', 'chowiluviennacl',
                                        'rowscalingviennacl', 'saviennacl',
                                        'bddc', 'kaczmarz', 'telescope']

        petsc_lu_direct_solvers = ['mumps', 'superlu_dist', 'umfpack', 'klu']

        petsc_cholesky_direct_solvers = ['mumps', 'cholmod']

        if self.preconditioner_type not in petsc_preconditioner_methods:
            self.preconditioner_type = 'jacobi'

        if ((self.solver_type in petsc_lu_direct_solvers) or
                (self.solver_type in petsc_cholesky_direct_solvers)):
            self.ksp = PETSc.KSP()
            self.ksp.create(PETSc.COMM_WORLD)

            self.ksp.getPC().setType(self.factorization_type)
            self.ksp.getPC().setFactorSolverPackage(self.solver_type)
            self.ksp.setType('preonly')

        elif self.solver_type in petsc_iterative_solvers:
            self.ksp = PETSc.KSP()
            self.ksp.create(PETSc.COMM_WORLD)

            self.ksp.getPC().setType(self.preconditioner_type)
            self.ksp.setType(self.solver_type)
        self.ksp.setTolerances(atol=self.atol,
                               rtol=self.rtol,
                               max_it=self.max_it)

    def _initialize_petsc_rhs_and_sol_vecs(self):
        r"""
        Initialize the solution vector (self.petsc_x), which is a dense
        matrix (1D vector) and defines the rhs vector (self.petsc_b) from
        the existing data.
        """

        # Get vector(s) compatible with the matrix,
        # i.e., with the same parallel layout.
        self.petsc_x, self.petsc_b = self.petsc_A.getVecs()

        #  Set the solution vector to zeros.
        self.petsc_x.set(0)

        # Define the petsc rhs vector from the numpy one.
        # If the rhs is defined by blocks, use this:
        PETSc.Vec.setValuesBlocked(self.petsc_b, [sp.arange(self.m)], self.b)
        # Otherwise, use:
        # PETSc.Vec.createWithArray(self.petsc_b,self.b)

    def petsc_solve(self, solver='cg',
                    preconditioner='jacobi',
                    factorization_type='lu',
                    atol=1e-06,
                    rtol=1e-06,
                    max_it=500):
        r"""
        This method solves the sparse linear system, converts the
        solution vector from a PETSc.Vec instance to a numpy array,
        and finally destroys all the petsc objects to free memory.

        Parameters
        ----------
        solver_type : string, optional
            Default is the iterative solver 'cg' based on the
            Conjugate Gradient method.
        preconditioner_type : string, optional
            Default is the 'jacobi' preconditioner, i.e., diagonal
            scaling preconditioning. The preconditioner is used with
            iterative solvers. When a direct solver is used, this
            parameter is ignored.
        factorization_type : string, optional
            The factorization type used with the direct solver.
            Default is 'lu'. This parameter is ignored when an
            iterative solver is used.

        Returns
        -------
        Returns a numpy array corresponding to the solution of
        the linear sparse system Ax = b.

        Examples
        --------
        >>> import scipy as sp
        >>> from openpnm_petsc_sparse_linear_solvers import \
        petsc_sparse_linear_solver as psls
        >>> ls = psls(alg.A, alg.b)
        >>> sol=psls.petsc_solve(ls)

        Notes
        -----
        Certain combinations of iterative solvers and precondioners
        or direct solvers and factorization types are not supported.
        The summary table of the different possibilities
        can be found here:
        https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html
        """
        self.solver_type = solver
        self.preconditioner_type = preconditioner
        self.factorization_type = factorization_type
        self.atol = atol
        self.rtol = rtol
        self.max_it = max_it

        self._initialize_petsc_sparse_coef_mat()
        self._create_petsc_sparse_linear_solver()
        self._initialize_petsc_rhs_and_sol_vecs()

        # PETSc
        self.ksp.setOperators(self.petsc_A)
        self.ksp.setFromOptions()

        self.ksp.solve(self.petsc_b, self.petsc_x)

        # Convert solution vector from PETSc.Vec instance
        # to a numpy array
        self.solution = PETSc.Vec.getArray(self.petsc_x)

        # Destroy petsc solver, coefficients matrix, rhs, and solution vectors
        PETSc.KSP.destroy(self.ksp)
        PETSc.Mat.destroy(self.petsc_A)
        PETSc.Vec.destroy(self.petsc_b)
        PETSc.Vec.destroy(self.petsc_x)

        return(self.solution)
