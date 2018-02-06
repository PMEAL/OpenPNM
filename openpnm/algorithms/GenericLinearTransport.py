"""
===============================================================================
module GenericLinearTransport: Class for solving linear transport processes
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.csgraph as spgr
from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
logger = logging.getLogger(__name__)


class GenericLinearTransport(GenericAlgorithm):
    r"""
    """

    def __init__(self, network, phase, settings={}, **kwargs):
        self.settings.update({'phase': phase.name,
                              'conductance': None,
                              'quantity': None,
                              'solver': 'spsolve',
                              'sources': [],
                              'tolerance': 0.001})
        self.settings.update(settings)
        super().__init__(Np=network.Np, Nt=network.Nt, network=network, **kwargs)
        # Initialize label 'all' in the object's own info dictionaries
        self['pore._id'] = network['pore._id']
        self['throat._id'] = network['throat._id']

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

    def set_source(self, source):
        self.settings['sources'].append(source.name)
        source.setup(algorithm=self)

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

        values = sp.array(bcvalues)
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
            self['pore.'+bctype+'_value'] = sp.nan
        self['pore.'+bctype+'_value'][pores] = values

    def remove_BC(self, pores=None):
        r"""
        Removes all boundary conditions assigned to the specified pores

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
            self['pore.dirichlet_value'][pores] = sp.nan
        if 'pore.neumann' in self.keys():
            self['pore.neumann'][pores] = False
            self['pore.neumann_value'][pores] = sp.nan

    def setup(self, conductance, quantity):
        r"""
        Specifies the necessary parameters

        Parameters
        ----------
        conductance : string
            The dictionary key containing the calculated pore-scale
            conductances.  For example, for StokesFlow this is
            'throat.hydraulic_conductance' by default.

        quantity : string
            The dictionary key where the values computed by this algorithm are
            stored.  For exaple, for StokesFLow this is 'pore.pressure' by
            default.

        """
        # Providing conductance values for the algorithm from the Physics name
        self.settings['conductance'] = self._parse_prop(conductance, 'throat')
        self.settings['quantity'] = self._parse_prop(quantity, 'pore')

        # Check health of conductance vector
        phase = self.simulation.phases[self.settings['phase']]
        if sp.any(sp.isnan(phase[self.settings['conductance']])):
            raise Exception('The provided throat conductance contains NaNs')

    def build_A(self):
        r"""
        """
        network = self.simulation.network
        phase = self.simulation.phases[self.settings['phase']]
        g = phase[self.settings['conductance']]
        am = network.create_adjacency_matrix(weights=-g, fmt='coo')
        A = spgr.laplacian(am)
        if 'pore.neumann' in self.keys():
            pass  # Do nothing to A, only b changes
        if 'pore.dirichlet' in self.keys():
            # Find all entries on rows associated with dirichlet pores
            P_bc = self.toindices(self['pore.dirichlet'])
            indrow = sp.in1d(A.row, P_bc)
            A.data[indrow] = 0  # Remove entries from A for all BC rows
            datadiag = A.diagonal()  # Add diagonal entries back into A
            datadiag[P_bc] = sp.ones_like(P_bc, dtype=float)
            A.setdiag(datadiag)
            A.eliminate_zeros()  # Remove 0 entries
        self.A = A
        return A

    def build_b(self):
        r"""
        """
        b = sp.zeros(shape=(self.Np, ), dtype=float)  # Create b matrix of 0's
        if 'pore.dirichlet' in self.keys():
            ind = self['pore.dirichlet']
            b[ind] = -self['pore.dirichlet_value'][ind]
        if 'pore.neumann' in self.keys():
            ind = self['pore.neumann']
            b[ind] = -self['pore.neumann_value'][ind]
        self.b = b
        return b

    def run(self, x=None):
        r"""
        Builds the A and b matrices, and calls the solver specified in the
        ``settings`` attribute.

        Parameters
        ----------
        x : ND-array
            Initial guess of unknown variable

        """
        self.build_A()
        self.build_b()
        if x is None:
            x = sp.zeros(shape=[self.Np, ], dtype=float)
        # Scan through all source objects registered on algorithm
        for item in self.settings['sources']:
            source = self.simulation[item]  # Obtain handle to source object
            source.apply()  # Apply source object to update A and b
        x_new = self.solve()
        self[self.settings['quantity']] = x_new
        res = sp.sum(sp.absolute(x**2 - x_new**2))
        if res < self.settings['tolerance']:
            print('Tolerance met, solution converged')
            return
        else:
            print('Tolerance not met: ', res)
            self.run(x=x_new)

    def solve(self, A=None, b=None):
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
        network = self.simulation.network
        phase = self.simulation.phases[self.settings['phase']]
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
            throat_group_num = sp.shape(t)[0]
        for i in sp.r_[0: throat_group_num]:
            if mode == 'group':
                throats = t
                P = pores
            elif mode == 'single':
                throats = t[i]
                P = pores[i]
            p1 = network.find_connected_pores(throats)[:, 0]
            p2 = network.find_connected_pores(throats)[:, 1]
            pores1 = sp.copy(p1)
            pores2 = sp.copy(p2)
            # Changes to pores1 and pores2 to make them as inner/outer pores
            pores1[~sp.in1d(p1, P)] = p2[~sp.in1d(p1, P)]
            pores2[~sp.in1d(p1, P)] = p1[~sp.in1d(p1, P)]
            X1 = quantity[pores1]
            X2 = quantity[pores2]
            g = conductance[throats]
            R.append(sp.sum(sp.multiply(g, (X2 - X1))))
        return sp.array(R, ndmin=1)

    def _calc_eff_prop(self):
        r"""
        Returns the main parameters for calculating the effective
        property in a linear transport equation.  It also checks for the
        proper boundary conditions, inlets and outlets.
        """
        network = self.simulation.network
        if self.settings['quantity'] not in self.keys():
            raise Exception('The algorit hm has not been run yet. Cannot ' +
                            'calculate effective property.')

        # Determine boundary conditions by analyzing algorithm object
        Ps = self.pores('pore.dirichlet')
        BCs = sp.unique(self['pore.dirichlet_value'][Ps])
        inlets = sp.where(self['pore.dirichlet_value'] == sp.amax(BCs))[0]
        outlets = sp.where(self['pore.dirichlet_value'] == sp.amin(BCs))[0]

        # Fetch area and length of domain
        A = network.domain_area(face=inlets)
        L = network.domain_length(face_1=inlets, face_2=outlets)
        flow = self.rate(pores=inlets)
        D = sp.sum(flow)*L/A/(BCs[0] - BCs[1])
        return D
