import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as spgr
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree
from openpnm.topotools import iscoplanar
from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
# check if petsc4py is available
import importlib
if (importlib.util.find_spec('petsc4py') is not None):
    from openpnm.utils.petsclinsolv import petscSparseLinearSolver as sls
logger = logging.getLogger(__name__)


class GenericTransport(GenericAlgorithm):
    r"""
    """

    def __init__(self, project=None, network=None, phase=None, settings={},
                 **kwargs):
        # Set some default settings
        self.settings.update({'phase': None,
                              'conductance': None,
                              'quantity': None,
                              'solver': 'spsolve'})
        self.setup(phase=phase, **settings)
        # If network given, get project, otherwise let parent class create it
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)
        # Create some instance attributes
        self._A = None
        self._pure_A = None
        self._b = None
        self._pure_b = None

    def setup(self, phase=None, **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Notes
        -----
        This generic version should be subclassed, and the arguments given
        suitable default names.
        """
        if phase:
            self.settings['phase'] = phase.name
        self.settings.update(kwargs)

    def set_value_BC(self, pores, values):
        r"""
        Apply constant value boundary conditons to the specified pore
        locations. These are sometimes referred to as Dirichlet conditions.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied

        values : scalar or array_like
            The value to of the boundary condition.  If a scalar is supplied
            it is assigne to all locations, and if a vector is applied it
            corresponds directy to the locations given in ``pores``.

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quentity'] = 'pore.pressure'``.
        """
        self._set_BC(pores=pores, bctype='value', bcvalues=values,
                     mode='merge')

    def set_rate_BC(self, pores, values):
        r"""
        Apply constant rate boundary conditons to the specified pore
        locations. This is similar to a Neumann boundary condition, but is
        slightly different since it's the conductance multiplied by the
        gradient, while Neumann conditions specify just the gradient.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied

        values : scalar or array_like
            The value to of the boundary condition.  If a scalar is supplied
            it is assigne to all locations, and if a vector is applied it
            corresponds directy to the locations given in ``pores``.

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quentity'] = 'pore.pressure'``.
        """
        self._set_BC(pores=pores, bctype='rate', bcvalues=values, mode='merge')

    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        Apply boundary conditions to specified pores

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied

        bctype : string
            Specifies the type or the name of boundary condition to apply. The
            types can be one one of the following:

            - *'value'* : Specify the value of the quantity in each location
            - *'rate'* : Specify the flow rate into each location

        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.  If
            a single value is given, it's assumed to apply to all locations.
            Different values can be applied to all pores in the form of an
            array of the same length as ``pores``.

        mode : string, optional
            Controls how the conditions are applied.  Options are:

            - *'merge'*: (Default) Adds supplied boundary conditions to already
            existing conditions.

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm. Use ``mode='remove'`` to
        clear existing BCs before applying new ones or ``mode='overwrite'``
        which removes all existing BC's before applying the new ones.

        Instead of using ``mode='remove'`` you can also set certain locations
        to NaN using ``mode='merge'``, which is equivalent to removing the BCs
        from those locations.

        """
        # Hijack the parse_mode function to verify bctype argument
        bctype = self._parse_mode(bctype, allowed=['value', 'rate'],
                                  single=True)
        mode = self._parse_mode(mode, allowed=['merge', 'overwrite', 'remove'],
                                single=True)
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of boundary values must match the ' +
                            'number of locations')

        # Store boundary values
        if ('pore.bc_'+bctype not in self.keys()) or (mode == 'overwrite'):
            self['pore.bc_'+bctype] = np.nan
        self['pore.bc_'+bctype][pores] = values

    def remove_BC(self, pores=None):
        r"""
        Removes all boundary conditions from the specified pores

        Parameters
        ----------
        pores : array_like
            The pores from which boundary conditions are to be removed.  If no
            pores are specified, then BCs are removed from all pores. No error
            is thrown if the provided pores do not have any BCs assigned.
        """
        if pores is None:
            pores = self.Ps
        if 'pore.bc_value' in self.keys():
            self['pore.bc_value'][pores] = np.nan
        if 'pore.rate' in self.keys():
            self['pore.bc_rate'][pores] = np.nan

    def _build_A(self, force=False):
        r"""
        Builds the coefficient matrix based on conductances between pores.
        The conductance to use is specified in the algorithm's ``settings``
        under ``conductance``.  In subclasses (e.g. ``FickianDiffusion``)
        this is set by default, though it can be overwritten.

        Parameters
        ----------
        force : Boolean (default is ``False)
            If set to ``True`` then the A matrix is built from new.  If
            ``False`` (the default), a cached version of A is returned.  The
            cached version is *clean* in the sense that no boundary conditions
            or sources terms have been added to it.
        """
        if force:
            self._pure_A = None
        if self._pure_A is None:
            network = self.project.network
            phase = self.project.phases()[self.settings['phase']]
            g = phase[self.settings['conductance']]
            am = network.create_adjacency_matrix(weights=g, fmt='coo')
            self._pure_A = spgr.laplacian(am)
        self.A = self._pure_A.copy()

    def _build_b(self, force=False):
        r"""
        Builds the RHS matrix, without applying any boundary conditions or
        source terms. This method is trivial an basically creates a column
        vector of 0's.

        Parameters
        ----------
        force : Boolean (default is ``False``)
            If set to ``True`` then the b matrix is built from new.  If
            ``False`` (the default), a cached version of b is returned.  The
            cached version is *clean* in the sense that no boundary conditions
            or sources terms have been added to it.
        """
        if force:
            self._pure_b = None
        if self._pure_b is None:
            b = np.zeros(shape=(self.Np, ), dtype=float)  # Create vector of 0s
            self._pure_b = b
        self.b = self._pure_b.copy()

    def _get_A(self):
        if self._A is None:
            self._build_A(force=True)
        return self._A

    def _set_A(self, A):
        self._A = A

    A = property(fget=_get_A, fset=_set_A)

    def _get_b(self):
        if self._b is None:
            self._build_b(force=True)
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
            f = np.amax(np.absolute(self.A.data))
            # Update b
            ind = np.isfinite(self['pore.bc_value'])
            self.b[ind] = f*self['pore.bc_value'][ind]
            # Update A
            # Find all entries on rows associated with value bc
            P_bc = self.toindices(np.isfinite(self['pore.bc_value']))
            indrow = np.in1d(self.A.row, P_bc)
            self.A.data[indrow] = 0  # Remove entries from A for all BC rows
            datadiag = self.A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = f*np.ones_like(P_bc, dtype=float)
            self.A.setdiag(datadiag)
            self.A.eliminate_zeros()  # Remove 0 entries

    def run(self):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        Returns
        -------
        Nothing is returned...the solution is stored on the objecxt under
        ``pore.quantity`` where *quantity* is specified in the ``settings``
        attribute.

        """
        print('―'*80)
        print('Running GenericTransport')
        self._run_generic()

    def _run_generic(self):
        self._apply_BCs()
        x_new = self._solve()
        self[self.settings['quantity']] = x_new

    def _solve(self, A=None, b=None):
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

        Notes
        -----
        The solver used here is specified in the ``settings`` attribute of the
        algorithm.

        """
        if A is None:
            A = self.A
            if A is None:
                raise Exception('The A matrix has not been built yet')
        if b is None:
            b = self.b
            if b is None:
                raise Exception('The b matrix has not been built yet')

        if self.settings['solver'] == 'petsc':
            # Check if petsc is available
            petsc = importlib.util.find_spec('petsc4py')
            if not petsc:
                raise Exception('petsc is not installed')
            if not self.settings['petsc_solver']:
                self.settings['petsc_solver'] = 'cg'
            if not self.settings['petsc_precond']:
                self.settings['petsc_precond'] = 'jacobi'
            # Define the petsc linear system converting the scipy objects
            ls = sls(A=A.tocsr(), b=b)
            x = sls.petsc_solve(ls, solver=self.settings['petsc_solver'],
                                preconditioner=self.settings['petsc_precond'])
            del(ls)  # Clean
        else:
            solver = getattr(sprs.linalg, self.settings['solver'])
            x = solver(A=A.tocsr(), b=b)
        return x

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
            Controls how to return the rate.  Options are:

            *'group'*: (default) Returns the cumulative rate of material
            moving into the given set of pores

            *'single'* : Calculates the rate for each pore individually

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

        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        g = phase[self.settings['conductance']]
        quantity = self[self.settings['quantity']]

        P12 = network['throat.conns']
        X12 = quantity[P12]
        f = (-1)**np.argsort(X12, axis=1)[:, 1]
        Dx = np.abs(np.diff(X12, axis=1).squeeze())
        Qt = -f*g*Dx

        if len(throats) and len(pores):
            raise Exception('Must specify either pores or throats, not both')
        elif len(throats):
            R = np.absolute(Qt[throats])
            if mode == 'group':
                R = np.sum(R)
        elif len(pores):
            Qp = np.zeros((self.Np, ))
            np.add.at(Qp, P12[:, 0], -Qt)
            np.add.at(Qp, P12[:, 1], Qt)
            R = Qp[pores]
            if mode == 'group':
                R = np.sum(R)
        return np.array(R, ndmin=1)

    def _calc_eff_prop(self):
        r"""
        Returns the main parameters for calculating the effective property
        in a linear transport equation.  It also checks for the proper
        boundary conditions, inlets and outlets.
        """
        if self.settings['quantity'] not in self.keys():
            raise Exception('The algorithm has not been run yet. Cannot ' +
                            'calculate effective property.')

        # Determine boundary conditions by analyzing algorithm object
        inlets, outlets = self._get_inlets_and_outlets()
        Ps = np.isfinite(self['pore.bc_value'])
        BCs = np.unique(self['pore.bc_value'][Ps])
        Dx = np.abs(np.diff(BCs))

        # Fetch area and length of domain
        A = self.domain_area
        L = self.domain_length
        flow = self.rate(pores=inlets)
        D = np.sum(flow)*L/A/Dx
        return D

    def _get_inlets_and_outlets(self):
        # Determine boundary conditions by analyzing algorithm object
        Ps = np.isfinite(self['pore.bc_value'])
        BCs = np.unique(self['pore.bc_value'][Ps])
        inlets = np.where(self['pore.bc_value'] == np.amax(BCs))[0]
        outlets = np.where(self['pore.bc_value'] == np.amin(BCs))[0]
        return (inlets, outlets)

    def _get_domain_area(self):
        if not hasattr(self, '_area'):
            logger.warning('Attempting to estimate inlet area...will be low')
            network = self.project.network
            Pin, Pout = self._get_inlets_and_outlets()
            inlets = network['pore.coords'][Pin]
            outlets = network['pore.coords'][Pout]
            if not iscoplanar(inlets):
                raise Exception('Detected inlet pores are not coplanar')
            if not iscoplanar(outlets):
                raise Exception('Detected outlet pores are not coplanar')
            Nin = np.ptp(inlets, axis=0) > 0
            if Nin.all():
                raise Exception('Detected inlets are not oriented along a ' +
                                'principle axis')
            Nout = np.ptp(outlets, axis=0) > 0
            if Nout.all():
                raise Exception('Detected outlets are not oriented along a ' +
                                'principle axis')
            hull_in = ConvexHull(points=inlets[:, Nin])
            hull_out = ConvexHull(points=outlets[:, Nout])
            if hull_in.volume != hull_out.volume:
                raise Exception('Inlet and outlet faces are different area')
            self._area = hull_in.volume  # In 2D volume=area, area=perimeter
        return self._area

    def _set_domain_area(self, area):
        self._area = area

    domain_area = property(fget=_get_domain_area, fset=_set_domain_area)

    def _get_domain_length(self):
        if not hasattr(self, '_length'):
            logger.warning('Attempting to estimate domain length... ' +
                           'could be low if boundary pores were not added')
            network = self.project.network
            Pin, Pout = self._get_inlets_and_outlets()
            inlets = network['pore.coords'][Pin]
            outlets = network['pore.coords'][Pout]
            if not iscoplanar(inlets):
                raise Exception('Detected inlet pores are not coplanar')
            if not iscoplanar(outlets):
                raise Exception('Detected inlet pores are not coplanar')
            tree = cKDTree(data=inlets)
            Ls = np.unique(np.around(tree.query(x=outlets)[0], decimals=5))
            if np.size(Ls) != 1:
                raise Exception('A unique value of length could not be found')
            self._length = Ls[0]
        return self._length

    def _set_domain_length(self, length):
        self._length = length

    domain_length = property(fget=_get_domain_length, fset=_set_domain_length)
