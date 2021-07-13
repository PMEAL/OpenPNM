import numpy as np
import openpnm as op
import scipy.sparse.linalg
import warnings
from numpy.linalg import norm
import scipy.sparse.csgraph as spgr
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree
from openpnm.topotools import iscoplanar, is_fully_connected
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import logging, Docorator, GenericSettings, prettify_logger_message
# Uncomment this line when we stop supporting Python 3.6
# from dataclasses import dataclass, field
# from typing import List

docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='GenericTransportSettings',
                     sections=['Parameters', 'Other Parameters'])
@docstr.dedent
# Uncomment this line when we stop supporting Python 3.6
# @dataclass
class GenericTransportSettings(GenericSettings):
    r"""
    Defines the settings for GenericTransport algorithms

    Parameters
    ----------
    phase : str
        The name of the phase on which the algorithm acts
    quantity : str
        The name of the physical quantity to be calculated
    conductance : str
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------
    solver_family : str (default = 'scipy')
        The solver package to use.  OpenPNM currently supports ``scipy``,
        ``pyamg`` and ``petsc`` (if you have it installed).
    solver_type : str
        The specific solver to use.  For instance, if ``solver_family`` is
        ``scipy`` then you can specify any of the iterative solvers such as
        ``cg`` or ``gmres``. [More info here]
        (https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html),
    solver_preconditioner : str (default = ``jacobi``)
        Used by the PETSc solver to specify which preconditioner to use.
    solver_tol : float (default = 1e-8)
        Used to control the accuracy to which the iterative solver aims to
        achieve before stopping. Can roughly be interpreted as the number of
        significant digits you want in your solution, i.e. 1e-8 -> 8
    solver_atol : float
        Absolute tolerance as a stopping criterion when using iterative
        solvers, defined as:
            atol = norm(Ax-b) @ x_final
        Don't specify this parameter unless you know what you're doing. The
        algorithm automatically calculates it based on ``solver_tol``.
    solver_rtol : float
        Relative tolerance as a stopping criterion when using iterative
        solvers, defined as the ratio of the residual computed using the
        accepted solution to that using the initial guess, i.e.:
            rtol = norm(Ax-b) @ x_final / norm(Ax-b) @ x0
        Don't specify this parameter unless you know what you're doing. The
        algorithm automatically calculates it based on ``solver_tol``.
    solver_max_iter : int (default = 5000)
        Maximum number of iterations allowed when using iterative solvers.
    variable_props : list
        List of pore/throat properties whose values might change during the
        algorithm and thus, need to be iterated for the solution to converge,
        e.g. "pore.diffusivity" if diffusivity is concentration-dependent.
    cache_A : bool
        If ``True``, A matrix is cached and reused rather than getting rebuilt.
    cache_b : bool
        If ``True``, b vector is cached and reused rather than getting rebuilt.

    """

    phase = None
    conductance = None
    quantity = None
    solver_family = 'pypardiso'
    solver_type = 'spsolve'
    solver_preconditioner = 'jacobi'
    solver_tol = 1e-8
    solver_atol = None
    solver_rtol = None
    solver_max_iter = 5000
    cache_A = True
    cache_b = True


@docstr.get_sections(base='GenericTransport', sections=['Parameters'])
@docstr.dedent
class GenericTransport(GenericAlgorithm):
    r"""
    This class implements steady-state linear transport calculations

    Parameters
    ----------
    %(GenericAlgorithm.parameters)s

    Notes
    -----

    The following table shows the methods that are accessible to the user
    for setting up the simulation.

    +---------------------+---------------------------------------------------+
    | Methods             | Description                                       |
    +=====================+===================================================+
    | ``set_value_BC``    | Applies constant value boundary conditions to the |
    |                     | specified pores                                   |
    +---------------------+---------------------------------------------------+
    | ``set_rate_BC``     | Applies constant rate boundary conditions to the  |
    |                     | specified pores                                   |
    +---------------------+---------------------------------------------------+
    | ``remove_BC``       | Removes all boundary conditions from the          |
    |                     | specified pores                                   |
    +---------------------+---------------------------------------------------+
    | ``rate``            | Calculates the total rate of transfer through the |
    |                     | given pores or throats                            |
    +---------------------+---------------------------------------------------+
    | ``setup``           | A shortcut for applying values in the ``settings``|
    |                     | attribute.                                        |
    +---------------------+---------------------------------------------------+
    | ``results``         | Returns the results of the calcualtion as a       |
    |                     | ``dict`` with the data stored under the 'quantity'|
    |                     | specified in the ``settings``                     |
    +---------------------+---------------------------------------------------+

    In addition to the above methods there are also the following attributes:

    +---------------------+---------------------------------------------------+
    | Attribute           | Description                                       |
    +=====================+===================================================+
    | ``A``               | Retrieves the coefficient matrix                  |
    +---------------------+---------------------------------------------------+
    | ``b``               | Retrieves the RHS matrix                          |
    +---------------------+---------------------------------------------------+

    This class contains quite a few hidden methods (preceeded by an
    underscore) that are called internally.  Since these are critical to the
    functioning of this algorithm they are worth outlining even though the
    user does not call them directly:

    +-----------------------+-------------------------------------------------+
    | Method or Attribute   | Description                                     |
    +=======================+=================================================+
    | ``_build_A``          | Builds the **A** matrix based on the            |
    |                       | 'conductance' specified in ``settings``         |
    +-----------------------+-------------------------------------------------+
    | ``_build_b``          | Builds the **b** matrix                         |
    +-----------------------+-------------------------------------------------+
    | ``_apply_BCs``        | Applies the given BCs by adjust the **A** and   |
    |                       | **b** matrices                                  |
    +-----------------------+-------------------------------------------------+
    | ``_calc_eff_prop``    | Finds the effective property (e.g. permeability |
    |                       | coefficient) based on the given BCs             |
    +-----------------------+-------------------------------------------------+
    | ``_solve``            | Runs the algorithm using the solver specified   |
    |                       | in the ``settings``                             |
    +-----------------------+-------------------------------------------------+
    | ``_get_domain_area``  | Attempts to estimate the area of the inlet pores|
    |                       | if not specified by user                        |
    +-----------------------+-------------------------------------------------+
    | ``_get_domain_length``| Attempts to estimate the length between the     |
    |                       | inlet and outlet faces if not specified by the  |
    |                       | user                                            |
    +-----------------------+-------------------------------------------------+

    """
    def __new__(cls, *args, **kwargs):
        instance = super(GenericTransport, cls).__new__(cls, *args, **kwargs)
        # Create some instance attributes
        instance._A = None
        instance._b = None
        instance._pure_A = None
        instance._pure_b = None
        return instance

    def __init__(self, project=None, network=None, phase=None, settings={},
                 **kwargs):
        # Apply default settings
        self.settings._update_settings_and_docs(GenericTransportSettings)
        # Overwrite any given in init
        self.settings.update(settings)
        # Assign phase if given during init
        self.setup(phase=phase)
        # If network given, get project, otherwise let parent class create it
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)
        self['pore.bc_rate'] = np.nan
        self['pore.bc_value'] = np.nan

    @docstr.get_sections(base='GenericTransport.setup', sections=['Parameters'])
    @docstr.dedent
    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        Customize algorithm settings, e.g. assign phase, set quantity to be
        solved, set conductance dict key, etc.

        Parameters
        ----------
        %(GenericTransportSettings.parameters)s

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        self.settings.update(**kwargs)

    @docstr.get_full_description(base='GenericTransport.reset')
    @docstr.get_sections(base='GenericTransport.reset', sections=['Parameters'])
    @docstr.dedent
    def reset(self, bcs=False, results=True):
        r"""
        Resets the algorithm to enable re-use.

        This allows the reuse of an algorithm inside a for-loop for parametric
        studies.  The default behavior means that only ``alg.reset()`` and
        ``alg.run()`` must be called inside a loop.  To reset the algorithm
        more completely requires overriding the default arguments.

        Parameters
        ----------
        results : boolean
            If ``True`` (default) all previously calculated values pertaining
            to results of the algorithm are removed.
        bcs : boolean (default = ``False``)
            If ``True`` all previous boundary conditions are removed.

        """
        self._pure_b = None
        self._b = None
        self._pure_A = None
        self._A = None
        if bcs:
            self['pore.bc_value'] = np.nan
            self['pore.bc_rate'] = np.nan
        if results:
            self.pop(self.settings['quantity'], None)

    @docstr.dedent
    def set_value_BC(self, pores, values, mode='merge'):
        r"""
        Apply constant value boundary conditons to the specified locations.

        These are sometimes referred to as Dirichlet conditions.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        values : scalar or array_like
            The value to apply in each pore.  If a scalar is supplied
            it is assigne to all locations, and if a vector is applied is
            must be the same size as the indices given in ``pores``.
        mode : string, optional
            Controls how the boundary conditions are applied.  Options are:

            'merge' - (Default) Adds supplied boundary conditions to already
            existing conditions, and also overwrites any existing values.
            If BCs of the complementary type already exist in the given
            locations, those values are kept.
            'overwrite' - Deletes all boundary conditions of the given type
            then adds the specified new ones (unless locations already have
            BCs of the other type).

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quantity'] = 'pore.pressure'``.

        """
        self._set_BC(pores=pores, bctype='value', bcvalues=values, mode=mode)

    def set_rate_BC(self, pores, rates=None, total_rate=None, mode='merge',
                    **kwargs):
        r"""
        Apply constant rate boundary conditons to the specified locations.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        rates : scalar or array_like, optional
            The rates to apply in each pore.  If a scalar is supplied that
            rate is assigned to all locations, and if a vector is supplied
            it must be the same size as the indices given in ``pores``.
        total_rate : float, optional
            The total rate supplied to all pores.  The rate supplied by this
            argument is divided evenly among all pores. A scalar must be
            supplied! Total_rate cannot be specified if rate is specified.
        mode : str, optional
            Controls how the boundary conditions are applied.  Options are:

            'merge' - (Default) Adds supplied boundary conditions to already
            existing conditions, and also overwrites any existing values.
            If BCs of the complementary type already exist in the given
            locations, these values are kept.
            'overwrite' - Deletes all boundary conditions of the given type
            then adds the specified new ones (unless locations already have
            BCs of the other type).

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quantity'] = 'pore.pressure'``.

        """
        # support 'values' keyword
        if 'values' in kwargs.keys():
            rates = kwargs.pop("values")
            warnings.warn("'values' has been deprecated, use 'rates' instead.",
                          DeprecationWarning)
        # handle total_rate feature
        if total_rate is not None:
            if not np.isscalar(total_rate):
                raise Exception('total_rate argument accepts scalar only!')
            if rates is not None:
                raise Exception('Cannot specify both arguments: rate and '
                                + 'total_rate')
            pores = self._parse_indices(pores)
            rates = total_rate/pores.size
        self._set_BC(pores=pores, bctype='rate', bcvalues=rates, mode=mode)

    @docstr.get_sections(base='GenericTransport._set_BC',
                         sections=['Parameters', 'Notes'])
    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        This private method is called by public facing BC methods, to apply
        boundary conditions to specified pores

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied
        bctype : string
            Specifies the type or the name of boundary condition to apply. The
            types can be one one of the following:

            'value' - Specify the value of the quantity in each location
            'rate' - Specify the flow rate into each location

        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.  If
            a single value is given, it's assumed to apply to all locations
            unless the 'total_rate' bc_type is supplied whereby a single value
            corresponds to a total rate to be divded evenly among all pores.
            Otherwise, different values can be applied to all pores in the form
            of an array of the same length as ``pores``.
        mode : string, optional
            Controls how the boundary conditions are applied.  Options are:

            'merge' - (Default) Adds supplied boundary conditions to already
            existing conditions, and also overwrites any existing values.
            If BCs of the complementary type already exist in the given
            locations, these values are kept.
            'overwrite' - Deletes all boundary conditions of the given type
            then adds the specified new ones (unless locations already have
            BCs of the other type).

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm. Use ``remove_BCs`` to
        clear existing BCs before applying new ones or ``mode='overwrite'``
        which removes all existing BC's before applying the new ones.

        """
        # Hijack the parse_mode function to verify bctype argument
        bctype = self._parse_mode(bctype, allowed=['value', 'rate'],
                                  single=True)
        othertype = list(set(['value', 'rate']).difference(set([bctype])))[0]
        mode = self._parse_mode(mode, allowed=['merge', 'overwrite'],
                                single=True)
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of boundary values must match the '
                            + 'number of locations')

        # Create boundary array if needed (though these are created on init)
        if 'pore.bc_' + bctype not in self.keys():
            self['pore.bc_' + bctype] = np.nan

        # Catch pores with existing BCs
        if mode == 'merge':  # remove offenders, and warn user
            existing_bcs = np.isfinite(self["pore.bc_" + othertype])
            inds = pores[existing_bcs[pores]]
        elif mode == 'overwrite':  # Remove existing BCs and write new ones
            self['pore.bc_' + bctype] = np.nan
            existing_bcs = np.isfinite(self["pore.bc_" + othertype])
            inds = pores[existing_bcs[pores]]
        # Now drop any pore indices which have BCs that should be kept
        if len(inds) > 0:
            msg = (r'Boundary conditions are already specified in the following given'
                   f' pores, so these will be skipped: {inds.__repr__()}')
            logger.warning(prettify_logger_message(msg))
            pores = np.setdiff1d(pores, inds)

        # Store boundary values
        self['pore.bc_' + bctype][pores] = values

    def remove_BC(self, pores=None, bctype='all'):
        r"""
        Removes boundary conditions from the specified pores

        Parameters
        ----------
        pores : array_like, optional
            The pores from which boundary conditions are to be removed.  If no
            pores are specified, then BCs are removed from all pores. No error
            is thrown if the provided pores do not have any BCs assigned.
        bctype : string, or list of strings
            Specifies which type of boundary condition to remove. Options are:

            -*'all'*: (default) Removes all boundary conditions
            -*'value'*: Removes only value conditions
            -*'rate'*: Removes only rate conditions

        """
        if isinstance(bctype, str):
            bctype = [bctype]
        if 'all' in bctype:
            bctype = ['value', 'rate']
        if pores is None:
            pores = self.Ps
        if ('pore.bc_value' in self.keys()) and ('value' in bctype):
            self['pore.bc_value'][pores] = np.nan
        if ('pore.bc_rate' in self.keys()) and ('rate' in bctype):
            self['pore.bc_rate'][pores] = np.nan

    def _build_A(self):
        r"""
        Builds the coefficient matrix based on conductances between pores.
        The conductance to use is specified in the algorithm's ``settings``
        under ``conductance``.  In subclasses (e.g. ``FickianDiffusion``)
        this is set by default, though it can be overwritten.
        """
        gvals = self.settings['conductance']
        if not gvals:
            raise Exception('conductance has not been defined on this algorithm')
        # Decide if caching of A and b is allowed
        # FIXME: this needs to be properly addressed (see issue #1548)
        try:
            if gvals in self._get_iterative_props():
                self.settings.update({"cache_A": False, "cache_b": False})
        except AttributeError:
            pass
        if not self.settings['cache_A']:
            self._pure_A = None
        if self._pure_A is None:
            network = self.project.network
            try:
                phase = self.project.phases()[self.settings['phase']]
            except KeyError:
                raise Exception('Phase has not been defined for algorithm')
            g = phase[gvals]
            am = network.create_adjacency_matrix(weights=g, fmt='coo')
            self._pure_A = spgr.laplacian(am).astype(float)
        self.A = self._pure_A.copy()

    def _build_b(self):
        r"""
        Builds the RHS matrix, without applying any boundary conditions or
        source terms. This method is trivial an basically creates a column
        vector of 0's.
        """
        cache_b = self.settings['cache_b']
        if not cache_b:
            self._pure_b = None
        if self._pure_b is None:
            b = np.zeros(shape=self.Np, dtype=float)  # Create vector of 0s
            self._pure_b = b
        self.b = self._pure_b.copy()

    def _get_A(self):
        if self._A is None:
            self._build_A()
        return self._A

    def _set_A(self, A):
        self._A = A

    A = property(fget=_get_A, fset=_set_A)

    def _get_b(self):
        if self._b is None:
            self._build_b()
        return self._b

    def _set_b(self, b):
        self._b = b

    b = property(fget=_get_b, fset=_set_b)

    def _apply_BCs(self):
        r"""
        Applies all the boundary conditions that have been specified, by
        adding values to the *A* and *b* matrices.
        """
        if 'pore.bc_rate' in self.keys():
            # Update b
            ind = np.isfinite(self['pore.bc_rate'])
            self.b[ind] = self['pore.bc_rate'][ind]
        if 'pore.bc_value' in self.keys():
            f = self.A.diagonal().mean()
            # Update b (impose bc values)
            ind = np.isfinite(self['pore.bc_value'])
            self.b[ind] = self['pore.bc_value'][ind] * f
            # Update b (substract quantities from b to keep A symmetric)
            x_BC = np.zeros_like(self.b)
            x_BC[ind] = self['pore.bc_value'][ind]
            self.b[~ind] -= (self.A * x_BC)[~ind]
            # Update A
            P_bc = self.toindices(ind)
            mask = np.isin(self.A.row, P_bc) | np.isin(self.A.col, P_bc)
            self.A.data[mask] = 0  # Remove entries from A for all BC rows/cols
            datadiag = self.A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = np.ones_like(P_bc, dtype=float) * f
            self.A.setdiag(datadiag)
            self.A.eliminate_zeros()  # Remove 0 entries

    def run(self, x0=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x0 : ND-array
            Initial guess of unknown variable

        Returns
        -------
        Nothing is returned...the solution is stored on the objecxt under
        ``pore.quantity`` where *quantity* is specified in the ``settings``
        attribute.

        """
        logger.info('―' * 80)
        logger.info('Running GenericTransport')
        self._validate_settings()
        # Check if A and b are well-defined
        self._validate_data_health()
        x0 = np.zeros_like(self.b) if x0 is None else x0
        self["pore.initial_guess"] = x0
        self._run_generic(x0)

    def _run_generic(self, x0):
        # (Re)build A,b in case phase/physics are updated and alg.run()
        # is to be called a second time
        self._build_A()
        self._build_b()
        self._apply_BCs()
        x_new = self._solve(x0=x0)
        quantity = self.settings['quantity']
        if not quantity:
            raise Exception('"quantity" has not been defined on this algorithm')
        self[quantity] = x_new

    def _solve(self, A=None, b=None, x0=None):
        r"""
        Sends the A and b matrices to the specified solver, and solves for *x*
        given the boundary conditions, and source terms based on the present
        value of *x*.  This method does NOT iterate to solve for non-linear
        source terms or march time steps.

        Parameters
        ----------
        A : sparse matrix
            The coefficient matrix in sparse format. If not specified, then
            it uses  the ``A`` matrix attached to the object.
        b : ND-array
            The RHS matrix in any format.  If not specified, then it uses
            the ``b`` matrix attached to the object.
        x0 : ND-array
            The initial guess for the solution of Ax = b

        Notes
        -----
        The solver used here is specified in the ``settings`` attribute of the
        algorithm.

        """
        x0 = np.zeros_like(self.b) if x0 is None else x0

        # Fetch A and b from self if not given, and throw error if not found
        A = self.A if A is None else A
        b = self.b if b is None else b
        if A is None or b is None:
            raise Exception('The A matrix or the b vector not yet built.')
        A = A.tocsr()

        # Check if A and b are STILL well-defined
        self._validate_data_health()

        # Check if A is symmetric
        if self.settings['solver_type'] == 'cg':
            is_sym = op.utils.is_symmetric(self.A)
            if not is_sym:
                raise Exception('CG solver only works on symmetric matrices.')

        # Fetch additional parameters for iterative solvers
        max_it = self.settings["solver_max_iter"]
        atol = self._get_atol()
        rtol = self._get_rtol(x0=x0)

        # Fetch solver object based on settings dict.
        solver = self._get_solver()
        x = solver(A, b, atol=atol, rtol=rtol, max_it=max_it, x0=x0)

        # Check solution convergence
        if not self._is_converged(x=x):
            raise Exception("Solver did not converge.")

        return x

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

    def _get_atol(self):
        r"""
        Fetches absolute tolerance for the solver if not ``None``, otherwise
        calculates it in a way that meets the given ``tol`` requirements.

        Notes
        -----
        ``atol`` is defined such to satisfy the following stopping criterion:
            ``norm(A*x-b)`` <= ``atol``

        """
        atol = self.settings["solver_atol"]
        if atol is None:
            tol = self.settings["solver_tol"]
            atol = norm(self.b) * tol
        return atol

    def _get_rtol(self, x0):
        r"""
        Fetches relative tolerance for the solver if not ``None``, otherwise
        calculates it in a way that meets the given ``tol`` requirements.

        ``rtol`` is defined based on the following formula:
            ``rtol = residual(@x_final) / residual(@x0)``
        """
        rtol = self.settings["solver_rtol"]
        if rtol is None:
            res0 = self._get_residual(x=x0)
            atol = self._get_atol()
            rtol = atol / res0
        return rtol

    def _get_residual(self, x=None):
        r"""
        Calculate solution residual based on the given ``x`` based on the
        following formula:
            ``res = norm(A*x - b)``
        """
        if x is None:
            quantity = self.settings['quantity']
            x = self[quantity]
        return norm(self.A * x - self.b)

    def _is_converged(self, x=None):
        r"""
        Check if solution has converged based on the following criterion:
            res <= max(norm(b) * tol, atol)
        """
        res = self._get_residual(x=x)
        # Verify that residual is finite (i.e. not inf/nan)
        if not np.isfinite(res):
            logger.error(f'Solution diverged: {res:.4e}')
            raise Exception(f"Solution diverged, undefined residual: {res:.4e}")
        # Check convergence
        tol = self.settings["solver_tol"]
        res_tol = norm(self.b) * tol
        flag_converged = True if res <= res_tol else False
        return flag_converged

    def _validate_settings(self):
        if self.settings['quantity'] is None:
            raise Exception('"quantity" has not been defined on this algorithm')
        if self.settings['conductance'] is None:
            raise Exception('"conductance" has not been defined on this algorithm')

    def _validate_geometry_health(self):
        h = self.project.check_geometry_health()
        issues = []
        for k, v in h.items():
            if len(v) > 0:
                issues.append(k)
        if len(issues) > 0:
            raise Exception(
                r"Found the following critical issues with your geometry(ies):"
                f" {', '.join(issues)}. Run network.project.check_geometry_health() for"
                r" more details.")

    def _validate_topology_health(self):
        Ps = (self['pore.bc_rate'] > 0) + (self['pore.bc_value'] > 0)
        if not is_fully_connected(network=self.network, pores_BC=Ps):
            raise Exception(
                "Your network is clustered. Run h = net.check_network_health() followed"
                " by op.topotools.trim(net, pores=h['trim_pores']) to make your network"
                " fully connected.")

    def _validate_data_health(self):
        r"""
        Check whether A and b are well-defined, i.e. doesn't contain nans.
        """
        import networkx as nx
        from pandas import unique

        # Short-circuit subsequent checks if data are healthy
        if np.isfinite(self.A.data).all() and np.isfinite(self.b).all():
            return True
        # Validate network topology health
        self._validate_topology_health()
        # Validate geometry health
        self._validate_geometry_health()

        # Fetch phase/geometries/physics
        prj = self.network.project
        phase = prj.find_phase(self)
        geometries = prj.geometries().values()
        physics = prj.physics().values()

        # Locate the root of NaNs
        unaccounted_nans = []
        for geom, phys in zip(geometries, physics):
            objs = [phase, geom, phys]
            # Generate global dependency graph
            dg = nx.compose_all([x.models.dependency_graph(deep=True) for x in objs])
            d = {}  # maps prop -> obj.name
            for obj in objs:
                for k, v in obj.check_data_health().items():
                    if "Has NaNs" in v:
                        # FIXME: The next line doesn't cover multi-level props
                        base_prop = ".".join(k.split(".")[:2])
                        if base_prop in dg.nodes:
                            d[base_prop] = obj.name
                        else:
                            unaccounted_nans.append(base_prop)
            # Generate dependency subgraph for props with NaNs
            dg_nans = nx.subgraph(dg, d.keys())
            # Find prop(s)/object(s) from which NaNs have propagated
            root_props = [n for n in d.keys() if not nx.ancestors(dg_nans, n)]
            root_objs = unique([d[x] for x in nx.topological_sort(dg_nans)])
            # Throw error with helpful info on how to resolve the issue
            if root_props:
                raise Exception(
                    r"Found NaNs in A matrix, possibly caused by NaNs in"
                    f" {', '.join(root_props)}. The issue might get resolved if you call"
                    r" regenerate_models on the following object(s):"
                    f" {', '.join(root_objs)}")

        # Raise Exception for unaccounted properties
        if unaccounted_nans:
            raise Exception(
                r"Found NaNs in A matrix, possibly caused by NaNs in"
                f" {', '.join(unaccounted_nans)}.")

        # Raise Exception otherwise if root cannot be found
        raise Exception(
            "Found NaNs in A matrix but couldn't locate the root object(s) that might"
            " have caused it. It's likely that disabling caching of A matrix via"
            " alg.settings['cache_A'] = False after instantiating the algorithm object"
            " fixes the problem.")

    def results(self):
        r"""
        Fetches the calculated quantity from the algorithm and returns it as
        an array.
        """
        quantity = self.settings['quantity']
        d = {quantity: self[quantity]}
        return d

    def rate(self, pores=[], throats=[], mode='group'):
        r"""
        Calculates the net rate of material moving into a given set of pores or
        throats

        Parameters
        ----------
        pores : array_like
            The pores for which the rate should be calculated
        throats : array_like
            The throats through which the rate should be calculated
        mode : string, optional
            Controls how to return the rate. Options are:
            - *'group'*: (default) Returns the cumulative rate of material
            moving into the given set of pores
            - *'single'* : Calculates the rate for each pore individually

        Returns
        -------
        If ``pores`` are specified, then the returned values indicate the
        net rate of material exiting the pore or pores.  Thus a positive
        rate indicates material is leaving the pores, and negative values
        mean material is entering.

        If ``throats`` are specified the rate is calculated in the direction of
        the gradient, thus is always positive.

        If ``mode`` is 'single' then the cumulative rate through the given
        pores (or throats) are returned as a vector, if ``mode`` is 'group'
        then the individual rates are summed and returned as a scalar.

        """
        pores = self._parse_indices(pores)
        throats = self._parse_indices(throats)

        if throats.size > 0 and pores.size > 0:
            raise Exception('Must specify either pores or throats, not both')
        if throats.size == pores.size == 0:
            raise Exception('Must specify either pores or throats')

        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        g = phase[self.settings['conductance']]
        quantity = self[self.settings['quantity']]

        P12 = network['throat.conns']
        X12 = quantity[P12]
        if g.size == self.Nt:
            g = np.tile(g, (2, 1)).T    # Make conductance a Nt by 2 matrix
        # The next line is critical for rates to be correct
        g = np.flip(g, axis=1)
        Qt = np.diff(g*X12, axis=1).squeeze()

        if throats.size:
            R = np.absolute(Qt[throats])
            if mode == 'group':
                R = np.sum(R)

        if pores.size:
            Qp = np.zeros((self.Np, ))
            np.add.at(Qp, P12[:, 0], -Qt)
            np.add.at(Qp, P12[:, 1], Qt)
            R = Qp[pores]
            if mode == 'group':
                R = np.sum(R)

        return np.array(R, ndmin=1)

    def set_solver(
            self,
            solver_family=None,
            solver_type=None,
            preconditioner=None,
            tol=None,
            atol=None,
            rtol=None,
            max_iter=None,
    ):
        r"""
        Set the solver to be used to solve the algorithm.

        The values of those fields that are not provided will be retrieved from
        algorithm settings dict.

        Parameters
        ----------
        solver_family : string, optional
            Solver family, could be "scipy", "petsc", and "pyamg".
        solver_type : string, optional
            Solver type, could be "spsolve", "cg", "gmres", etc.
        preconditioner : string, optional
            Preconditioner for iterative solvers. The default is "jacobi".
        tol : float, optional
            Tolerance for iterative solvers, loosely related to number of
            significant digits in data.
        atol : float, optional
            Absolute tolerance for iterative solvers, such that
            norm(Ax-b) <= atol holds.
        rtol : float, optional
            Relative tolerance for iterative solvers, loosely related to how
            many orders of magnitude reduction in residual is desired, compared
            to its value at initial guess.
        max_iter : int, optional
            Maximum number of iterations

        Returns
        -------
        None

        """
        settings = self.settings
        # Preserve pre-set values, if any
        if solver_family is None:
            solver_family = settings["solver_family"]
        if solver_type is None:
            solver_type = settings["solver_type"]
        if preconditioner is None:
            preconditioner = settings["solver_preconditioner"]
        if tol is None:
            tol = settings["solver_tol"]
        if atol is None:
            atol = settings["solver_atol"]
        if rtol is None:
            rtol = settings["solver_rtol"]
        if max_iter is None:
            max_iter = settings["solver_max_iter"]
        # Update settings on algorithm object
        self.settings.update(
            {
                "solver_family": solver_family,
                "solver_type": solver_type,
                "solver_preconditioner": preconditioner,
                "solver_tol": tol,
                "solver_atol": atol,
                "solver_rtol": rtol,
                "solver_max_iter": max_iter
            }
        )

    def _calc_eff_prop(self, inlets=None, outlets=None,
                       domain_area=None, domain_length=None):
        r"""
        Calculate the effective transport through the network

        Parameters
        ----------
        inlets : array_like
            The pores where the inlet boundary conditions were applied.  If
            not given an attempt is made to infer them from the algorithm.
        outlets : array_like
            The pores where the outlet boundary conditions were applied.  If
            not given an attempt is made to infer them from the algorithm.
        domain_area : scalar
            The area of the inlet and/or outlet face (which shold match)
        domain_length : scalar
            The length of the domain between the inlet and outlet faces

        Returns
        -------
        The effective transport property through the network

        """
        if self.settings['quantity'] not in self.keys():
            raise Exception('The algorithm has not been run yet. Cannot '
                            + 'calculate effective property.')

        Ps = np.isfinite(self['pore.bc_value'])
        BCs = np.unique(self['pore.bc_value'][Ps])
        Dx = np.abs(np.diff(BCs))
        if inlets is None:
            inlets = self._get_inlets()
        flow = self.rate(pores=inlets)
        # Fetch area and length of domain
        if domain_area is None:
            domain_area = self._get_domain_area(inlets=inlets,
                                                outlets=outlets)
        if domain_length is None:
            domain_length = self._get_domain_length(inlets=inlets,
                                                    outlets=outlets)
        D = np.sum(flow)*domain_length/domain_area/Dx
        return D

    def _get_inlets(self):
        # Determine boundary conditions by analyzing algorithm object
        Ps = np.isfinite(self['pore.bc_value'])
        BCs = np.unique(self['pore.bc_value'][Ps])
        inlets = np.where(self['pore.bc_value'] == np.amax(BCs))[0]
        return inlets

    def _get_outlets(self):
        # Determine boundary conditions by analyzing algorithm object
        Ps = np.isfinite(self['pore.bc_value'])
        BCs = np.unique(self['pore.bc_value'][Ps])
        outlets = np.where(self['pore.bc_value'] == np.amin(BCs))[0]
        return outlets

    def _get_domain_area(self, inlets=None, outlets=None):
        logger.warning('Attempting to estimate inlet area...will be low')
        network = self.project.network
        # Abort if network is not 3D
        if np.sum(np.ptp(network['pore.coords'], axis=0) == 0) > 0:
            raise Exception('The network is not 3D, specify area manually')
        if inlets is None:
            inlets = self._get_inlets()
        if outlets is None:
            outlets = self._get_outlets()
        inlets = network['pore.coords'][inlets]
        outlets = network['pore.coords'][outlets]
        if not iscoplanar(inlets):
            logger.error('Detected inlet pores are not coplanar')
        if not iscoplanar(outlets):
            logger.error('Detected outlet pores are not coplanar')
        Nin = np.ptp(inlets, axis=0) > 0
        if Nin.all():
            logger.warning('Detected inlets are not oriented along a '
                           + 'principle axis')
        Nout = np.ptp(outlets, axis=0) > 0
        if Nout.all():
            logger.warning('Detected outlets are not oriented along a '
                           + 'principle axis')
        hull_in = ConvexHull(points=inlets[:, Nin])
        hull_out = ConvexHull(points=outlets[:, Nout])
        if hull_in.volume != hull_out.volume:
            logger.error('Inlet and outlet faces are different area')
        area = hull_in.volume  # In 2D volume=area, area=perimeter
        return area

    def _get_domain_length(self, inlets=None, outlets=None):
        logger.warning('Attempting to estimate domain length...'
                       + 'could be low if boundary pores were not added')
        network = self.project.network
        if inlets is None:
            inlets = self._get_inlets()
        if outlets is None:
            outlets = self._get_outlets()
        inlets = network['pore.coords'][inlets]
        outlets = network['pore.coords'][outlets]
        if not iscoplanar(inlets):
            logger.error('Detected inlet pores are not coplanar')
        if not iscoplanar(outlets):
            logger.error('Detected inlet pores are not coplanar')
        tree = cKDTree(data=inlets)
        Ls = np.unique(np.float64(tree.query(x=outlets)[0]))
        if not np.allclose(Ls, Ls[0]):
            logger.error('A unique value of length could not be found')
        length = Ls[0]
        return length
