import logging
import numpy as np
import scipy.sparse.csgraph as spgr
from openpnm.topotools import is_fully_connected
from openpnm.algorithms import Algorithm
from openpnm.utils import Docorator, TypedSet, Workspace
from openpnm.utils import check_data_health
from openpnm import solvers
from ._solution import SteadyStateSolution, SolutionContainer


__all__ = ['Transport']


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


@docstr.get_sections(base='TransportSettings', sections=['Parameters'])
@docstr.dedent
class TransportSettings:
    """
    Defines the settings for Transport algorithms

    Parameters
    ----------
    %(AlgorithmSettings.parameters)s
    quantity : str
        The name of the physical quantity to be solved for (i.e.
        'pore.concentration')
    conductance : str
        The name of the pore-scale transport conductance values (i.e
        'throat.diffusive_conductance')
    cache : bool
        If ``True``, the ``A`` matrix is cached and rather than getting
        rebuilt.
    variable_props : list of strings
        This list (actually a set) indicates which properties are variable
        and should be updated by the algorithm on each iteration. Note that
        any properties which already depend on ``'quantity'`` will
        automatically be updated.

    """
    phase = ''
    quantity = ''
    conductance = ''
    cache = True
    variable_props = TypedSet()


@docstr.get_sections(base='Transport', sections=['Parameters'])
@docstr.dedent
class Transport(Algorithm):
    """
    This class implements steady-state linear transport calculations.

    Parameters
    ----------
    %(Algorithm.parameters)s

    """

    def __init__(self, phase, name='trans_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(TransportSettings())
        self.settings['phase'] = phase.name
        self['pore.bc.rate'] = np.nan
        self['pore.bc.value'] = np.nan
        self._A = None
        self._b = None
        self._pure_A = None
        self._pure_b = None
        self.soln = {}

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            if key == 'pore.source':
                return {}
            else:
                raise KeyError(key)

    @property
    def x(self):
        """Shortcut to the solution currently stored on the algorithm."""
        return self[self.settings['quantity']]

    @x.setter
    def x(self, value):
        self[self.settings['quantity']] = value

    @docstr.dedent
    def _build_A(self):
        """
        Builds the coefficient matrix based on throat conductance values.

        Notes
        -----
        The conductance to use is specified in stored in the algorithm's
        settings under ``alg.settings['conductance']``.

        """
        gvals = self.settings['conductance']
        if gvals in self.iterative_props:
            self.settings.cache = False
        if not self.settings['cache']:
            self._pure_A = None
        if self._pure_A is None:
            phase = self.project[self.settings.phase]
            g = phase[gvals]
            am = self.network.create_adjacency_matrix(weights=g, fmt='coo')
            self._pure_A = spgr.laplacian(am).astype(float)
        self.A = self._pure_A.copy()

    def _build_b(self):
        """Initializes the RHS vector, b, with zeros."""
        b = np.zeros(self.Np, dtype=float)
        self._pure_b = b
        self.b = self._pure_b.copy()

    @property
    def A(self):
        """The coefficient matrix, A (in Ax = b)"""
        if self._A is None:
            self._build_A()
        return self._A

    @A.setter
    def A(self, value):
        self._A = value

    @property
    def b(self):
        """The right-hand-side (RHS) vector, b (in Ax = b)"""
        if self._b is None:
            self._build_b()
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    def _apply_BCs(self):
        """Applies specified boundary conditions by modifying A and b."""
        if 'pore.bc.rate' in self.keys():
            # Update b
            ind = np.isfinite(self['pore.bc.rate'])
            self.b[ind] = self['pore.bc.rate'][ind]
        if 'pore.bc.value' in self.keys():
            f = self.A.diagonal().mean()
            # Update b (impose bc values)
            ind = np.isfinite(self['pore.bc.value'])
            self.b[ind] = self['pore.bc.value'][ind] * f
            # Update b (subtract quantities from b to keep A symmetric)
            x_BC = np.zeros_like(self.b)
            x_BC[ind] = self['pore.bc.value'][ind]
            self.b[~ind] -= (self.A * x_BC)[~ind]
            # Update A
            P_bc = self.to_indices(ind)
            mask = np.isin(self.A.row, P_bc) | np.isin(self.A.col, P_bc)
            # Remove entries from A for all BC rows/cols
            self.A.data[mask] = 0
            # Add diagonal entries back into A
            datadiag = self.A.diagonal()
            datadiag[P_bc] = np.ones_like(P_bc, dtype=float) * f
            self.A.setdiag(datadiag)
            self.A.eliminate_zeros()

    def run(self, solver=None, x0=None, verbose=False):
        """
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        This method stores the solution in the algorithm's ``soln``
        attribute as a ``SolutionContainer`` object. The solution itself
        is stored in the ``x`` attribute of the algorithm as a NumPy array.

        Parameters
        ----------
        x0 : ndarray
            Initial guess of unknown variable

        Returns
        -------
        None

        """
        logger.info('Running Transport')
        if solver is None:
            solver = getattr(solvers, ws.settings.default_solver)()
        # Perform pre-solve validations
        self._validate_settings()
        self._validate_topology_health()
        self._validate_linear_system()
        # Write x0 to algorithm (needed by _update_iterative_props)
        self.x = x0 = np.zeros_like(self.b) if x0 is None else x0.copy()
        self["pore.initial_guess"] = x0
        self._validate_x0()
        # Initialize the solution object
        self.soln = SolutionContainer()
        self.soln[self.settings['quantity']] = SteadyStateSolution(x0)
        self.soln.is_converged = False
        # Build A and b, then solve the system of equations
        self._update_A_and_b()
        self._run_special(solver=solver, x0=x0, verbose=verbose)

    def _run_special(self, solver, x0, w=1.0, verbose=None):
        # Make sure A and b are 'still' well-defined
        self._validate_linear_system()
        # Solve and apply under-relaxation
        x_new, exit_code = solver.solve(A=self.A, b=self.b, x0=x0)
        self.x = w * x_new + (1 - w) * self.x
        # Update A and b using the recent solution otherwise, for iterative
        # algorithms, residual will be incorrectly calculated ~0, since A & b
        # are outdated
        self._update_A_and_b()
        # Update SteadyStateSolution object on algorithm
        self.soln[self.settings['quantity']][:] = self.x
        self.soln.is_converged = not bool(exit_code)

    def _update_A_and_b(self):
        """Builds A and b, and applies specified boundary conditions."""
        self._build_A()
        self._build_b()
        self._apply_BCs()

    def _validate_x0(self):
        """Ensures x0 doesn't contain any nans/infs."""
        x0 = self["pore.initial_guess"]
        if not np.isfinite(x0).all():
            raise Exception("x0 contains inf/nan values")

    def _validate_settings(self):
        if self.settings['quantity'] is None:
            raise Exception("'quantity' hasn't been defined on this algorithm")
        if self.settings['conductance'] is None:
            raise Exception("'conductance' hasn't been defined on this algorithm")
        if self.settings['phase'] is None:
            raise Exception("'phase' hasn't been defined on this algorithm")

    def _validate_topology_health(self):
        """
        Ensures the network is not clustered, and if it is, they're at
        least connected to a boundary condition pore.
        """
        Ps = ~np.isnan(self['pore.bc.rate']) + ~np.isnan(self['pore.bc.value'])
        if not is_fully_connected(network=self.network, pores_BC=Ps):
            msg = ("Your network is clustered, making Ax = b ill-conditioned")
            raise Exception(msg)

    def _validate_linear_system(self):
        """Ensures the linear system Ax = b doesn't contain any nans/infs."""
        if np.isfinite(self.A.data).all() and np.isfinite(self.b).all():
            return
        raise Exception("A or b contains inf/nan values")

    def rate(self, pores=[], throats=[], mode='group'):
        """
        Calculates the net rate of material moving into a given set of
        pores or throats

        Parameters
        ----------
        pores : array_like
            The pores for which the rate should be calculated
        throats : array_like
            The throats through which the rate should be calculated
        mode : str, optional
            Controls how to return the rate. The default value is 'group'.
            Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'group'      Returns the cumulative rate of material
            'single'     Calculates the rate for each pore individually
            ===========  =====================================================

        Returns
        -------
        If ``pores`` are specified, then the returned values indicate the
        net rate of material exiting the pore or pores.  Thus a positive
        rate indicates material is leaving the pores, and negative values
        mean material is entering.

        If ``throats`` are specified the rate is calculated in the
        direction of the gradient, thus is always positive.

        If ``mode`` is 'single' then the cumulative rate through the given
        pores (or throats) are returned as a vector, if ``mode`` is
        'group' then the individual rates are summed and returned as a
        scalar.

        """
        pores = self._parse_indices(pores)
        throats = self._parse_indices(throats)

        if throats.size > 0 and pores.size > 0:
            raise Exception('Must specify either pores or throats, not both')
        if (throats.size == 0) and (pores.size == 0):
            raise Exception('Must specify either pores or throats')

        network = self.project.network
        phase = self.project[self.settings['phase']]
        g = phase[self.settings['conductance']]

        P12 = network['throat.conns']
        X12 = self.x[P12]
        if g.size == self.Nt:
            g = np.tile(g, (2, 1)).T    # Make conductance an Nt by 2 matrix
        # The next line is critical for rates to be correct
        # We could also do "g.T.flatten()" or "g.flatten('F')"
        g = np.flip(g, axis=1)
        Qt = np.diff(g*X12, axis=1).ravel()

        if throats.size:
            R = np.absolute(Qt[throats])
            if mode == 'group':
                R = np.sum(R)
        elif pores.size:
            Qp = np.zeros((self.Np, ))
            np.add.at(Qp, P12[:, 0], -Qt)
            np.add.at(Qp, P12[:, 1], Qt)
            R = Qp[pores]
            if mode == 'group':
                R = np.sum(R)

        return np.array(R, ndmin=1)

    def clear_value_BCs(self):
        """Clears all value BCs."""
        self.set_BC(pores=None, bctype='value', mode='remove')

    def clear_rate_BCs(self):
        """Clears all rate BCs"""
        self.set_BC(pores=None, bctype='rate', mode='remove')

    def set_value_BC(self, pores=None, values=[], mode='add'):
        """
        Applies constant value boundary conditons to the specified pores.

        These are sometimes referred to as Dirichlet conditions.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        values : float or array_like
            The value to apply in each pore. If a scalar is supplied
            it is assigne to all locations, and if a vector is applied is
            must be the same size as the indices given in ``pores``.
        mode : str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. For definition of various modes, see the
            docstring for ``set_BC``.

        """
        self.set_BC(pores=pores, bctype='value', bcvalues=values, mode=mode)

    def set_rate_BC(self, pores=None, rates=[], mode='add'):
        """
        Apply constant rate boundary conditons to the specified locations.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        rates : float or array_like, optional
            The rates to apply in each pore. If a scalar is supplied that
            rate is assigned to all locations, and if a vector is supplied
            it must be the same size as the indices given in ``pores``.
        mode : str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. For definition of various modes, see the
            docstring for ``set_BC``.
        """
        self.set_BC(pores=pores, bctype='rate', bcvalues=rates, mode=mode)
