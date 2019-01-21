import importlib
import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as spgr
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree
from decimal import Decimal as dc
from openpnm.topotools import iscoplanar
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import logging
logger = logging.getLogger(__name__)

# Set some default settings
def_set = {'phase': None,
           'conductance': None,
           'quantity': None,
           'solver_family': 'scipy',
           'solver_type': 'spsolve',
           'solver_preconditioner': 'jacobi',
           'solver_atol': 1e-6,
           'solver_rtol': 1e-6,
           'solver_maxiter': 5000,
           'gui': {'setup':        {'quantity': '',
                                    'conductance': ''},
                   'set_rate_BC':  {'pores': None,
                                    'values': None},
                   'set_value_BC': {'pores': None,
                                    'values': None},
                   'remove_BC':    {'pores': None}
                   }
           }


class GenericTransport(GenericAlgorithm):
    r"""
    This class implements steady-state linear transport calculations

    Parameters
    ----------
    network : OpenPNM Network object
        The Network with which this algorithm is associated

    project : OpenPNM Project object, optional
        A Project can be specified instead of ``network``

    Notes
    -----

    The following table shows the methods that are accessible to the user
    for settig up the simulation.

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

    def __init__(self, project=None, network=None, phase=None, settings={},
                 **kwargs):
        # Apply default settings
        self.settings.update(def_set)
        # Overwrite any given in init
        self.settings.update(settings)
        # Assign phase if given during init
        self.setup(phase=phase)
        # If network given, get project, otherwise let parent class create it
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)
        # Create some instance attributes
        self._A = None
        self._pure_A = None
        self._b = None
        self._pure_b = None
        self['pore.bc_rate'] = np.nan
        self['pore.bc_value'] = np.nan

    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run.

        quantity : string
            The name of the physical quantity to be calculated.

        conductance : string
            The name of the pore-scale transport conductance values.  These
            are typically calculated by a model attached to a *Physics* object
            associated with the given *Phase*.

        solver : string
            To use the default scipy solver, set this value to `spsolve` or
            `umfpack`.  To use an iterative solver or a non-scipy solver,
            additional arguments are required as described next.

        solver_family : string
            The solver package to use.  OpenPNM currently supports ``scipy``,
            ``pyamg`` and ``petsc`` (if you have it installed).  The default is
            ``scipy``.

        solver_type : string
            The specific solver to use.  For instance, if ``solver_family`` is
            ``scipy`` then you can specify any of the iterative solvers such as
            ``cg`` or ``gmres``.  [More info here]
            (https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html)

        solver_preconditioner : string
            This is used by the PETSc solver to specify which preconditioner
            to use.  The default is ``jacobi``.

        solver_atol : scalar
            Used to control the accuracy to which the iterative solver aims.
            The default is 1e-6.

        solver_rtol : scalar
            Used by PETSc as an additional tolerance control.  The default is
            1e-6.

        solver_maxiter : scalar
            Limits the number of iterations to attempt before quiting when
            aiming for the specified tolerance. The default is 5000.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        self.settings.update(**kwargs)

    def set_value_BC(self, pores, values):
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

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quentity'] = 'pore.pressure'``.
        """
        self._set_BC(pores=pores, bctype='value', bcvalues=values,
                     mode='merge')

    def set_rate_BC(self, pores, values):
        r"""
        Apply constant rate boundary conditons to the specified locations.

        This is similar to a Neumann boundary condition, but is
        slightly different since it's the conductance multiplied by the
        gradient, while Neumann conditions specify just the gradient.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied

        values : scalar or array_like
            The values of rate to apply in each pore.  If a scalar is supplied
            it is assigned to all locations, and if a vector is applied it
            must be the same size as the indices given in ``pores``.

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quentity'] = 'pore.pressure'``.
        """
        self._set_BC(pores=pores, bctype='rate', bcvalues=values, mode='merge')

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

            - ``'value'``: Specify the value of the quantity in each location
            - ``'rate'``: Specify the flow rate into each location

        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.  If
            a single value is given, it's assumed to apply to all locations.
            Different values can be applied to all pores in the form of an
            array of the same length as ``pores``.

        mode : string, optional
            Controls how the boundary conditions are applied.  Options are:

            - ``'merge'``: (Default) Adds supplied boundary conditions to
            already existing conditions.
            - ``'overwrite'``: Deletes all boundary condition on object then
            adds the given ones

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
            Specifies which type of boundary condition to remove.  Options are:

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

    def _build_A(self, force=False):
        r"""
        Builds the coefficient matrix based on conductances between pores.
        The conductance to use is specified in the algorithm's ``settings``
        under ``conductance``.  In subclasses (e.g. ``FickianDiffusion``)
        this is set by default, though it can be overwritten.

        Parameters
        ----------
        force : Boolean (default is ``False``)
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
            f = np.abs(self.A.data).mean()
            # Update b (impose bc values)
            ind = np.isfinite(self['pore.bc_value'])
            self.b[ind] = self['pore.bc_value'][ind] * f
            # Update b (substract quantities from b to keep A symmetric)
            x_BC = np.zeros(self.b.shape)
            x_BC[ind] = self['pore.bc_value'][ind]
            self.b[~ind] -= (self.A.tocsr() * x_BC)[~ind]
            # Update A
            P_bc = self.toindices(ind)
            indrow = np.isin(self.A.row, P_bc)
            indcol = np.isin(self.A.col, P_bc)
            self.A.data[indrow] = 0  # Remove entries from A for all BC rows
            self.A.data[indcol] = 0  # Remove entries from A for all BC cols
            datadiag = self.A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = np.ones_like(P_bc, dtype=np.float64) * f
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
        logger.info('―'*80)
        logger.info('Running GenericTransport')
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
        # Fetch A and b from self if not given, and throw error if they've not
        # been calculated
        if A is None:
            A = self.A
            if A is None:
                raise Exception('The A matrix has not been built yet')
        if b is None:
            b = self.b
            if b is None:
                raise Exception('The b matrix has not been built yet')
        A = A.tocsr()

        # Default behavior -> use Scipy's default solver (spsolve)
        if self.settings['solver'] == 'pyamg':
            self.settings['solver_family'] = 'pyamg'
        if self.settings['solver'] == 'petsc':
            self.settings['solver_family'] = 'petsc'

        # Set tolerance for iterative solvers
        rtol = self.settings['solver_rtol']
        min_A = np.abs(A.data).min()
        min_b = np.abs(b).min() or 1e100
        atol = min(min_A, min_b) * rtol

        # SciPy
        if self.settings['solver_family'] == 'scipy':
            if importlib.util.find_spec('scikit-umfpack'):
                A.indices = A.indices.astype(np.int64)
                A.indptr = A.indptr.astype(np.int64)
            iterative = ['bicg', 'bicgstab', 'cg', 'cgs', 'gmres', 'lgmres',
                         'minres', 'gcrotmk', 'qmr']
            solver = getattr(sprs.linalg, self.settings['solver_type'])
            if self.settings['solver_type'] in iterative:
                x, exit_code = solver(A=A, b=b, atol=atol, tol=rtol,
                                      maxiter=self.settings['solver_maxiter'])
                if exit_code > 0:
                    raise Exception('SciPy solver did not converge! ' +
                                    'Exit code: ' + str(exit_code))
            else:
                x = solver(A=A, b=b)
            return x

        # PETSc
        if self.settings['solver_family'] == 'petsc':
            # Check if petsc is available
            if importlib.util.find_spec('petsc4py'):
                from openpnm.utils.petsc import PETScSparseLinearSolver as SLS
            else:
                raise Exception('PETSc is not installed.')
            # Define the petsc linear system converting the scipy objects
            ls = SLS(A=A, b=b)
            sets = self.settings
            sets = {k: v for k, v in sets.items() if k.startswith('solver_')}
            sets = {k.split('solver_')[1]: v for k, v in sets.items()}
            ls.settings.update(sets)
            x = SLS.solve(ls)
            del(ls)
            return x

        # PyAMG
        if self.settings['solver_family'] == 'pyamg':
            if importlib.util.find_spec('pyamg'):
                import pyamg
            else:
                raise Exception('pyamg is not installed.')
            ml = pyamg.ruge_stuben_solver(A)
            x = ml.solve(b=b, tol=1e-6)
            return x

    def results(self, times='all', t_precision=12, **kwargs):
        r"""
        Fetches the calculated quantity from the algorithm and returns it as
        an array.

        Parameters
        ----------
        times : scalar or list
            Time steps to be returned. The default value is 'all' which results
            in returning all time steps. If a scalar is given, only the
            corresponding time step is returned. If a range is given
            (e.g., 'range(0, 1, 1e-3)'), time steps in this range are returned.

        t_precision : integer
            The time precision (number of decimal places). Default value is 12.

        Notes
        -----
        The keyword steps is interpreted in the same way as times.
        """
        if 'steps' in kwargs.keys():
            times = kwargs['steps']
        t_pre = t_precision
        quantity = self.settings['quantity']
        q = [k for k in list(self.keys()) if quantity in k]
        if times == 'all':
            t = q
        elif type(times) in [float, int]:
            n = int(-dc(str(round(times, t_pre))).as_tuple().exponent *
                    (round(times, t_pre) != int(times)))
            t_str = (str(int(round(times, t_pre)*10**n)) +
                     ('e-'+str(n))*(n != 0))
            t = [k for k in q if t_str == k.split('@')[-1]]
        elif 'range' in times:
            t = times.replace(' ', '')
            t = t[6:-1]
            t = t.split(',')
            out = np.arange(float(t[0]), float(t[1]), float(t[2]))
            out = np.append(out, float(t[1]))
            out = np.unique(out)
            out = np.around(out, decimals=t_pre)
            t = []
            for i in out:
                n = int(-dc(str(round(i, t_pre))).as_tuple().exponent *
                        (round(i, t_pre) != int(i)))
                j = (str(int(round(i, t_pre)*10**n))+('e-'+str(n))*(n != 0))
                t_str = [k for k in q if j == k.split('@')[-1]]
                t += (t_str)
        d = {k: self[k] for k in t}
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
            raise Exception('The algorithm has not been run yet. Cannot ' +
                            'calculate effective property.')

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
            logger.warning('Detected inlets are not oriented along a ' +
                           'principle axis')
        Nout = np.ptp(outlets, axis=0) > 0
        if Nout.all():
            logger.warning('Detected outlets are not oriented along a ' +
                           'principle axis')
        hull_in = ConvexHull(points=inlets[:, Nin])
        hull_out = ConvexHull(points=outlets[:, Nout])
        if hull_in.volume != hull_out.volume:
            logger.error('Inlet and outlet faces are different area')
        area = hull_in.volume  # In 2D volume=area, area=perimeter
        return area

    def _get_domain_length(self, inlets=None, outlets=None):
        logger.warning('Attempting to estimate domain length... ' +
                       'could be low if boundary pores were not added')
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
        Ls = np.unique(np.around(tree.query(x=outlets)[0], decimals=5))
        if np.size(Ls) != 1:
            logger.error('A unique value of length could not be found')
        length = Ls[0]
        return length
