import logging
import numpy as np
import scipy.sparse.csgraph as spgr
from openpnm.topotools import is_fully_connected
from openpnm.algorithms import GenericAlgorithm
from openpnm.algorithms import BCsMixin
from openpnm.utils import Docorator, SettingsAttr, TypedSet, Workspace
from openpnm.utils import check_data_health
from openpnm import solvers
from ._solution import SteadyStateSolution, SolutionContainer


__all__ = ['GenericTransport', 'GenericTransportSettings']


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


@docstr.get_sections(base='GenericTransportSettings', sections=['Parameters'])
@docstr.dedent
class GenericTransportSettings:
    r"""
    Defines the settings for GenericTransport algorithms

    Parameters
    ----------
    %(GenericAlgorithmSettings.parameters)s
    quantity : str
        The name of the physical quantity to be calculated
    conductance : str
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.
    cache : bool
        If ``True``, A matrix is cached and rather than getting rebuilt.

    """
    prefix = 'transport'
    phase = ''
    quantity = ''
    conductance = ''
    cache = True
    variable_props = TypedSet()


@docstr.get_sections(base='GenericTransport', sections=['Parameters'])
@docstr.dedent
class GenericTransport(GenericAlgorithm, BCsMixin):
    r"""
    This class implements steady-state linear transport calculations.

    Parameters
    ----------
    %(GenericAlgorithm.parameters)s

    """

    def __new__(cls, *args, **kwargs):
        instance = super(GenericTransport, cls).__new__(cls, *args, **kwargs)
        # Create some instance attributes
        instance._A = None
        instance._b = None
        instance._pure_A = None
        instance._pure_b = None
        return instance

    def __init__(self, phase, settings=None, **kwargs):
        self.settings = SettingsAttr(GenericTransportSettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        self.settings['phase'] = phase.name
        self['pore.bc_rate'] = np.nan
        self['pore.bc_value'] = np.nan

    @property
    def x(self):
        """Shortcut to the solution currently stored on the algorithm."""
        return self[self.settings['quantity']]

    @x.setter
    def x(self, value):
        self[self.settings['quantity']] = value

    @docstr.dedent
    def _build_A(self):
        r"""
        Builds the coefficient matrix based on throat conductance values.

        Notes
        -----
        The conductance to use is specified in stored in the algorithm's
        settings under ``alg.settings['conductance']``.

        """
        gvals = self.settings['conductance']
        # FIXME: this needs to be properly addressed (see issue #1548)
        try:
            if gvals in self._get_iterative_props():
                self.settings.cache = False
        except AttributeError:
            pass
        if not self.settings['cache']:
            self._pure_A = None
        if self._pure_A is None:
            phase = self.project[self.settings.phase]
            g = phase[gvals]
            am = self.network.create_adjacency_matrix(weights=g, fmt='coo')
            self._pure_A = spgr.laplacian(am).astype(float)
        self.A = self._pure_A.copy()

    def _build_b(self):
        r"""
        Builds the RHS vector, without applying any boundary conditions or
        source terms. This method is trivial as it just creates a column
        vector of 0's.
        """
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
            # Update b (subtract quantities from b to keep A symmetric)
            x_BC = np.zeros_like(self.b)
            x_BC[ind] = self['pore.bc_value'][ind]
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

    def run(self, solver=None, x0=None, verbose=True):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x0 : ndarray
            Initial guess of unknown variable

        Notes
        -----
        This method doesn't return anything. The solution is stored on
        the object under ``pore.<quantity>`` where *<quantity>* is specified
        in the ``settings`` attribute.

        """
        logger.info('Running GenericTransport')
        if solver is None:
            solver = getattr(solvers, ws.settings.default_solver)()
        # Perform pre-solve validations
        self._validate_settings()
        self._validate_data_health()
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
        return self.soln

    def _run_special(self, solver, x0, w=1.0, verbose=None):
        # Make sure A,b are STILL well-defined
        self._validate_data_health()
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
        r"""
        Builds/updates A, b based on the recent solution on the algorithm
        object.
        """
        self._build_A()
        self._build_b()
        self._apply_BCs()

    def _validate_x0(self):
        """
        Ensures x0 doesn't contain any nans/infs.
        """
        x0 = self["pore.initial_guess"]
        if not np.isfinite(x0).all():
            raise Exception("x0 contains Inf/NaN values")

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
        Ps = ~np.isnan(self['pore.bc_rate']) + ~np.isnan(self['pore.bc_value'])
        if not is_fully_connected(network=self.network, pores_BC=Ps):
            msg = ("Your network is clustered. Run h = net.check_network_"
                   "health() followed by op.topotools.trim(net, pores=h['disconnected_pores']) to make your network fully connected.")
            raise Exception(msg)

    def _validate_conductance_model_health(self):
        """
        Ensures all throats have a conductance model assigned.
        """
        conductance = self.settings.conductance
        g = self.project.phases(self.settings.phase)[conductance]
        Ts_nan = self.Ts[np.isnan(g)]
        Ts_with_model = []
        for obj in self.project:
            if conductance in obj.keys():
                Ts_with_model.extend(obj.throats(to_global=True))
        if not np.in1d(Ts_nan, Ts_with_model).all():
            msg = ("Found NaNs in A matrix, possibly because some throats"
                   f" don't have conductance model assigned: {conductance}")
            raise Exception(msg)

    def _validate_data_health(self):
        """
        Checks whether A and b are well-defined, i.e. doesn't contain nans.
        """
        import networkx as nx
        from pandas import unique

        # Validate network topology health
        self._validate_topology_health()
        # Short-circuit subsequent checks if data are healthy
        if np.isfinite(self.A.data).all() and np.isfinite(self.b).all():
            return True

        # Fetch phase/geometries/physics
        prj = self.project
        phase = prj.phases(self.settings.phase)

        # Locate the root of NaNs
        unaccounted_nans = []
        objs = [phase]
        # Generate global dependency graph
        dg = nx.compose_all([x.models.dependency_graph(deep=True) for x in objs])
        d = {}  # maps prop -> obj.name
        for obj in objs:
            for k, v in check_data_health(obj).items():
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
            msg = ("Found NaNs in A matrix, possibly caused by NaNs in"
                   f" {', '.join(root_props)}. The issue might get resolved"
                   " if you call `regenerate_models` on the following"
                   f" object(s): {', '.join(root_objs)}")
            raise Exception(msg)

        # Raise Exception for throats without an assigned conductance model
        self._validate_conductance_model_health()

        # Raise Exception for unaccounted properties
        if unaccounted_nans:
            msg = ("Found NaNs in A matrix, possibly caused by NaNs in"
                   f" {', '.join(unaccounted_nans)}.")
            raise Exception(msg)

        # Raise Exception otherwise if root cannot be found
        msg = ("Found NaNs in A matrix but couldn't locate the root cause."
               " It's likely that disabling caching of A matrix via"
               " `alg.settings['cache'] = False` after instantiating the"
               " algorithm object fixes the problem.")
        raise Exception(msg)

    def rate(self, pores=[], throats=[], mode='group'):
        r"""
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
        g = np.flip(g, axis=1)
        Qt = np.diff(g*X12, axis=1).squeeze()

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

    def set_variable_props(self, variable_props, mode='add'):
        r"""
        This method is useful for setting variable_props to the settings
        dictionary of the target object. Variable_props and their dependent
        properties get updated iteratively.

        Parameters
        ----------
        variable_props : str, or List(str)
            A single string or list of strings to be added as variable_props
        mode : str, optional
            Controls how the variable_props are applied. The default value is
            'add'. Options are:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'add'        Adds supplied variable_props to already existing list
                         (if any), and prevents duplicates
            'overwrite'  Deletes all exisitng variable_props and then adds
                         the specified new ones
            ===========  =====================================================

        """
        # If single string, make it a list
        if isinstance(variable_props, str):
            variable_props = [variable_props]
        # Handle mode
        mode = self._parse_mode(mode, allowed=['add', 'overwrite'],
                                single=True)
        if mode == 'overwrite':
            self.settings['variable_props'].clear()
        # parse each propname and append to variable_props in settings
        for variable_prop in variable_props:
            variable_prop = self._parse_prop(variable_prop, 'pore')
            self.settings['variable_props'].add(variable_prop)
