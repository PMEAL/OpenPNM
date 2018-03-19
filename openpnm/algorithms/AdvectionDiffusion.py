import numpy as np
import scipy.sparse as sprs
import scipy.sparse.csgraph as spgr
from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
logger = logging.getLogger(__name__)


class AdvectionDiffusion(GenericAlgorithm):
    r"""
    A subclass of GenericAlgorithm to simulate advection diffusion.

    """

    def __init__(self, project=None, network=None, phase=None, settings={},
                 **kwargs):
        # Set some default settings
        self.settings.update({'phase': None,
                       'quantity': 'pore.mole_fraction',
                       'diffusive_conductance': 'throat.diffusive_conductance',
                       'hydraulic_conductance': 'throat.hydraulic_conductance',
                       'pressure': 'pore.pressure',
                       'solver': 'spsolve'})
        self.settings.update(settings)
        if phase is not None:
            self.settings['phase'] = phase.name

        if network is not None:
            project = network.project

        super().__init__(project=project, **kwargs)

    def set_dirchlet_BC(self, pores, values):
        r"""
        """
        self.set_BC(pores=pores, bctype='dirichlet', bcvalues=values,
                    mode='merge')

    def set_neumann_BC(self, pores, values):
        r"""
        """
        self.set_BC(pores=pores, bctype='neumann', bcvalues=values,
                    mode='merge')

    def set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        Apply boundary conditions to specified pores

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied

        bctype : string
            Specifies the type or the name of boundary condition to apply. The
            types can be one one of the following:

            - *'dirichlet'* : Specify the quantity in each location
            - *'neumann'* : Specify the flow rate into each location

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
        specified location in just one algorithm. Use ``mode='remove'`` to
        clear existing BCs before applying new ones or ``mode='overwrite'``
        which removes all existing BC's before applying the new ones.

        Instead of using ``mode='remove'`` you can also set certain locations
        to NaN using ``mode='merge'``, which is equivalent to removing the BCs
        from those locations.

        """
        # Hijack the parse_mode function to verify bctype argument
        bctype = self._parse_mode(bctype, allowed=['dirichlet', 'neumann',
                                                   'neumann_group'],
                                  single=True)
        mode = self._parse_mode(mode, allowed=['merge', 'overwrite', 'remove'],
                                single=True)
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of boundary values must match the ' +
                            'number of locations')

        # Label pores where a boundary condition will be applied
        if ('pore.'+bctype not in self.keys()) or (mode == 'overwrite'):
            self['pore.'+bctype] = False
        self['pore.'+bctype][pores] = True

        # Store boundary values
        if ('pore.'+bctype+'_value' not in self.keys()) or \
           (mode == 'overwrite'):
            self['pore.'+bctype+'_value'] = np.nan
        self['pore.'+bctype+'_value'][pores] = values

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
        if 'pore.dirichlet' in self.keys():
            self['pore.dirichlet'][pores] = False
            self['pore.dirichlet_value'][pores] = np.nan
        if 'pore.neumann' in self.keys():
            self['pore.neumann'][pores] = False
            self['pore.neumann_value'][pores] = np.nan

    def build_A(self):
        r"""
        """
        network = self.project.network
        phase = self.project.phases()[self.settings['phase']]
        HC = phase[self.settings['hydraulic_conductance']]
        p = phase[self.settings['pressure']]
        D = phase[self.settings['diffusive_conductance']]
        nt = network.find_neighbor_throats(pores=phase['pore.all'], \
             flatten=False, mode='not_intersection') # Pore neighbor throats
        A = np.zeros((network.Np, network.Np)) # Initialize A matrix
        for i in range (network.Np):
            q = HC[nt[i]]*(p[network.find_neighbor_pores(i)]-p[i]) # Flow rate
            qP = np.where(q>0, q, 0) # Throat positive flow rates
            qN = np.where(q<0, q, 0) # Throat negative flow rates
            A[i,i] = np.sum( qN - D[nt[i]] ) # Diagonal
            j1 = network['throat.conns'][nt[i]] # Find off diag cells to fill
            j2 = np.reshape(j1,np.size(j1))
            j = j2[j2!=i]
            A[i,j] = -( -qP - D[nt[i]] ) # Off diagonal
        A = sprs.coo_matrix(A)
        self.A = A
        return A

    def build_b(self):
        r"""
        """
        b = np.zeros(shape=(self.Np, ), dtype=float)  # Create b matrix of 0's
        self.b = b
        return b

    def apply_BCs(self):
        if 'pore.neumann' in self.keys():
            # Update b
            ind = self['pore.neumann']
            self.b[ind] = self['pore.neumann_value'][ind]

        if 'pore.dirichlet' in self.keys():
            # Update b
            ind = self['pore.dirichlet']
            self.b[ind] = self['pore.dirichlet_value'][ind]
            # Update A
            # Find all entries on rows associated with dirichlet pores
            P_bc = self.toindices(self['pore.dirichlet'])
            indrow = np.in1d(self.A.row, P_bc)
            self.A.data[indrow] = 0  # Remove entries from A for all BC rows
            datadiag = self.A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = np.ones_like(P_bc, dtype=float)
            self.A.setdiag(datadiag)
            self.A.eliminate_zeros()  # Remove 0 entries

    def setup(self):
        r"""
        """
        self.build_A()
        self.build_b()
        self.apply_BCs()

    def run(self):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        """
        print('―'*80)
        print('Running GenericTransport')
        self.setup()
        self._run_generic()

    def _run_generic(self):
        x_new = self._solve()
        self[self.settings['quantity']] = x_new

    def _solve(self, A=None, b=None):
        r"""
        Sends the A and b matrices to the specified solver.

        Parameters
        ----------
        A : sparse matrix
            The coefficient matrix in sparse format

        b : ND-array
            The RHS matrix in any format

        Notes
        -----
        The solver used here is specified in the ``settings`` attribute of the
        algorithm.

        """
        if A is None:
            A = self.A
        if b is None:
            b = self.b
        solver = getattr(sprs.linalg, self.settings['solver'])
        x = solver(A=A.tocsr(), b=b)
        return x

    def results(self):
        r"""
        """
        quantity = self.settings['quantity']
        d = {quantity: self[quantity]}
        return d

    def rate(self, pores=None, mode='group'):
        r"""
        Calculates the net rate of material moving into a given set of pores.

        Parameters
        ----------
        pores : array_like
            The pores for which the rate should be calculated

        mode : string, optional
            Controls how to return the rate.  Options are:

            *'group'*: (default) Teturns the cumulative rate of material
            moving into the given set of pores

            *'single'* : Calculates the rate for each pore individually

        Notes
        -----
        A negative rate indicates material moving into the pore or pores, such
        as material being consumed.
        """
        network = self.project.network
        phase = self.project.phases[self.settings['phase']]
        conductance = phase[self.settings['conductance']]
        quantity = self[self.settings['quantity']]
        pores = self._parse_indices(pores)
        R = []
        if mode == 'group':
            t = network.find_neighbor_throats(pores, flatten=True,
                                              mode='not_intersection')
            throat_group_num = 1
        elif mode == 'single':
            t = network.find_neighbor_throats(pores, flatten=False,
                                              mode='not_intersection')
            throat_group_num = np.shape(t)[0]
        for i in np.r_[0: throat_group_num]:
            if mode == 'group':
                throats = t
                P = pores
            elif mode == 'single':
                throats = t[i]
                P = pores[i]
            p1 = network.find_connected_pores(throats)[:, 0]
            p2 = network.find_connected_pores(throats)[:, 1]
            pores1 = np.copy(p1)
            pores2 = np.copy(p2)
            # Changes to pores1 and pores2 to make them as inner/outer pores
            pores1[~np.in1d(p1, P)] = p2[~np.in1d(p1, P)]
            pores2[~np.in1d(p1, P)] = p1[~np.in1d(p1, P)]
            X1 = quantity[pores1]
            X2 = quantity[pores2]
            g = conductance[throats]
            R.append(np.sum(np.multiply(g, (X2 - X1))))
        return np.array(R, ndmin=1)

    def _calc_eff_prop(self):
        r"""
        Returns the main parameters for calculating the effective
        property in a linear transport equation.  It also checks for the
        proper boundary conditions, inlets and outlets.
        """
        network = self.project.network
        if self.settings['quantity'] not in self.keys():
            raise Exception('The algorithm has not been run yet. Cannot ' +
                            'calculate effective property.')

        # Determine boundary conditions by analyzing algorithm object
        Ps = self.pores('pore.dirichlet')
        BCs = np.unique(self['pore.dirichlet_value'][Ps])
        inlets = np.where(self['pore.dirichlet_value'] == np.amax(BCs))[0]
        outlets = np.where(self['pore.dirichlet_value'] == np.amin(BCs))[0]

        # Fetch area and length of domain
        A = network.domain_area(face=inlets)
        L = network.domain_length(face_1=inlets, face_2=outlets)
        flow = self.rate(pores=inlets)
        D = np.sum(flow)*L/A/(BCs[0] - BCs[1])
        return D
