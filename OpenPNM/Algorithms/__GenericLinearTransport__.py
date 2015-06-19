# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericLinearTransport__: Class for solving linear transport processes
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as sprslin
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Phases import GenericPhase
import OpenPNM.Utilities.vertexops as vo
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class GenericLinearTransport(GenericAlgorithm):
    r"""
    This class provides essential methods for building and solving matrices
    in a transport process.  It is inherited by FickianDiffusion,
    FourierConduction, StokesFlow and OhmicConduction.
    """
    def __init__(self, phase=None, **kwargs):
        super().__init__(**kwargs)
        if phase is None:
            self._phase = GenericPhase()
            self._phases.append(self._phase)
        else:
            self._phase = phase  # Register phase with self
            if sp.size(phase) != 1:
                self._phases = phase
            else:
                self._phases.append(phase)
        for comp in self._phases:
            if comp.Np != self.Np:
                raise Exception(comp.name + ' has different Np size than the' +
                                ' algorithm ' + self.name)

    def setup(self, conductance, quantity, super_pore_conductance):
        r"""
        This setup provides the initial data for the solver from the provided
        properties.
        It also creates the matrices A and b.
        """
        # Assigning super_pore conductance for Neumann_group BC
        if super_pore_conductance is None:
            self.super_pore_conductance = []
        else:
            self.super_pore_conductance = super_pore_conductance
        # Providing conductance values for the algorithm from the Physics name
        if sp.size(self._phase) == 1:
            self._conductance = 'throat.'+conductance.split('.')[-1]
            self._quantity = 'pore.' + self._phase.name + '_' + \
                             quantity.split('.')[-1]
            # Check health of conductance vector
            if self._phase.check_data_health(props=self._conductance).health:
                self['throat.conductance'] = self._phase[self._conductance]
            else:
                raise Exception('The provided throat conductance has problems')
        else:
            raise Exception('The linear solver accepts just one phase.')
        # Checking for the linear terms to be added to the coeff diagonal/RHS
        diag_added_data = sp.zeros(self.Np)
        RHS_added_data = sp.zeros(self.Np)
        for label in self.labels():
            if 'pore.source_' in label:
                source_name = 'pore.' + \
                              (label.split('.')[-1]).replace('source_', '')
                matching_physics = [phys for phys in self._phase._physics
                                    if source_name in phys.models.keys()]
                for phys in matching_physics:
                    x = phys.models[source_name]['x']
                    if x != '' and type(x) == str:
                        if x.split('.')[-1] != quantity.split('.')[-1]:
                            raise Exception('The quantity(pore.' +
                                            x.split('.')[-1] +
                                            '), provided by source term(' +
                                            source_name + '), is different ' +
                                            'from the main quantity(pore.' +
                                            quantity.split('.')[-1] + ') in ' +
                                            self.name + ' algorithm.')
                source_name = label.replace('pore.source_', '')
                if 'pore.source_linear_s1_' + source_name in self.props():
                    prop1 = 'pore.source_linear_s1_' + source_name
                    pores = ~sp.isnan(self[prop1])
                    diag_added_data[pores] = diag_added_data[pores] + \
                        self[prop1][pores]
                    prop2 = 'pore.source_linear_s2_' + source_name
                    pores = ~sp.isnan(self[prop2])
                    RHS_added_data[pores] = RHS_added_data[pores] + \
                        self[prop2][pores]
        # Creating A and b based on the conductance values and new linear terms
        logger.info('Creating Coefficient matrix for the algorithm')
        d = diag_added_data
        self.A = self._build_coefficient_matrix(modified_diag_pores=self.Ps,
                                                diag_added_data=d)
        logger.info('Creating RHS matrix for the algorithm')
        self.b = self._build_RHS_matrix(modified_RHS_pores=self.Ps,
                                        RHS_added_data=-RHS_added_data)

    def set_source_term(self, source_name=None, pores=None, x0=None, tol=None,
                        maxiter=None, mode='merge'):
        r"""
        Apply source terms to specified pores

        Parameters
        ----------
        source_name : string
          Specifies the name of source term from a Physics object to apply.
        pores : array_like
          The pores where the boundary conditions should be applied
        x0 : array_like, optional
          By sending guess values for the quantity, the method calculates the
          source terms and stores them in the algorithm
        tol : float, optional
          Tolerance for the iterative method. (if maxiter>0)
        maxiter : integer, optional
          Maximum number of iterations for this source term. Iteration will
          stop after maxiter steps.
        mode : string, optional
          Controls how the source terms should be applied.
          Options are:
                - 'merge': Inserts specified values, leaving existing values \
                  elsewhere.
                - 'overwrite': Inserts specified values, clearing all other \
                  values.
                - 'remove': Removes boundary conditions from specified \
                  locations.
                - 'update': Allows to insert specified values to new \
                  locations, updating existing ones.

        Notes
        -----
        Difference between 'merge' and 'update' modes: in the merge, a new
        value cannot be applied to a pore with existing one, but in the
        'update' it is possible.

        """

        if mode not in ['merge', 'overwrite', 'remove', 'update']:
            raise Exception('The mode (' + mode + ') cannot be applied to ' +
                            'the set_source_term!')
        if pores is not None:
            pores = sp.array(pores, ndmin=1)
        # Checking for existance of source_name
        if source_name is not None:
            s_group = sp.array(source_name, ndmin=1)
            for source_name in s_group:
                source_name = 'pore.' + source_name.split('.')[-1]
                prop = source_name.split('.')[-1]
                try:
                    self._phase[source_name]
                except KeyError:
                    Exception('The attached phase in the algorithm ' +
                              self.name + ', does not have the source ' +
                              'property ' + source_name + ' in its physics!')
                except ValueError:
                    pass

                if mode == 'remove':
                    s_mode = ['linear', 'nonlinear']
                    if source_name is None:
                        if pores is not None:
                            if pores is 'all':
                                for item in self.labels():
                                    if 'pore.source_' in item:
                                        temp_item = item.split('.')[-1]
                                        prop = temp_item.replace('source_', '')
                                        del self['pore.source_' +
                                                 prop]
                                        for s in s_mode:
                                            try:
                                                del self['pore.source_' +
                                                         s + '_s1_' + prop]
                                            except KeyError:
                                                pass
                                            try:
                                                del self['pore.source_' +
                                                         s + '_s2_' + prop]
                                            except KeyError:
                                                pass
                            else:
                                for item in self.labels():
                                    if 'pore.source_' in item:
                                        temp_item = item.split('.')[-1]
                                        prop = temp_item.replace('source_', '')
                                        self['pore.source_' +
                                             prop][pores] = False
                                        for s in s_mode:
                                            try:
                                                self['pore.source_' +
                                                     s + '_s1_' +
                                                     prop][pores] = sp.nan
                                            except KeyError:
                                                pass
                                            try:
                                                self['pore.source_' +
                                                     s + '_s2_' +
                                                     prop][pores] = sp.nan
                                            except KeyError:
                                                pass
                        else:
                            raise Exception('No pores/source_name are ' +
                                            'sent to the set_term_method!')
                    else:
                        if pores is None:
                            try:
                                del self['pore.source_' + prop]
                            except KeyError:
                                pass
                            for s in s_mode:
                                try:
                                    del self['pore.source_' +
                                             s + '_s1_' + prop]
                                except KeyError:
                                    pass
                                try:
                                    del self['pore.source_' +
                                             s + '_s2_' + prop]
                                except KeyError:
                                    pass
                        else:
                            try:
                                self['pore.source_' + prop][pores] = False
                            except KeyError:
                                pass
                            for s in s_mode:
                                try:
                                    self['pore.source_' +
                                         s + '_s1_' +
                                         prop][pores] = sp.nan
                                except KeyError:
                                    pass
                                try:
                                    self['pore.source_' +
                                         s + '_s2_' + prop][pores] = sp.nan
                                except KeyError:
                                    pass

                else:
                    # Handle tol, x0 and maxiter for the Picard algorithm
                    if 'pore.source_tol' not in self.props():
                        self['pore.source_tol'] = sp.ones((self.Np,),
                                                          dtype=float) * sp.nan
                    if 'pore.source_maxiter' not in self.props():
                        maxiter_arr = sp.ones((self.Np,), dtype=float) * sp.nan
                        self['pore.source_maxiter'] = maxiter_arr

                    if x0 is None:
                        x0 = 0
                    self._guess = x0
                    # Check value of maxiter
                    if maxiter is None:
                        maxiter = int(100)
                        source_mode = 'nonlinear'
                    else:
                        try:
                            maxiter = int(maxiter)
                        except (ValueError, TypeError):
                            raise Exception('input for maxiter cannot be ' +
                                            'converted to integer!')
                        if maxiter > 0:
                            source_mode = 'nonlinear'
                        elif maxiter == 0:
                            source_mode = 'linear'
                    # Check value of tol
                    if tol is None:
                        tol = 1e-5
                    else:
                        try:
                            tol = float(tol)
                        except (ValueError, TypeError):
                            raise Exception('input for tol cannot be ' +
                                            'converted to float!')

                    if ('pore.source_' + prop not in self.labels() or
                            mode == 'overwrite'):
                        self['pore.source_' + prop] = sp.zeros((self.Np,),
                                                               dtype=bool)
                        temp_arr = sp.ones((self.Np,), dtype=float) * sp.nan
                        self['pore.source_' + source_mode +
                             '_s1_' + prop] = temp_arr
                        self['pore.source_' + source_mode +
                             '_s2_' + prop] = temp_arr

                    # Setting the source term for all the modes except 'remove'
                    matching_physics = [phys for phys in self._phase._physics
                                        if source_name in phys.models.keys()]
                    for phys in matching_physics:
                        x = phys.models[source_name]['x']
                        return_rate = phys.models[source_name]['return_rate']
                        regen_mode = phys.models[source_name]['regen_mode']
                        phys.models[source_name]['x'] = x0
                        phys.models[source_name]['return_rate'] = False
                        phys.models[source_name]['regen_mode'] = 'normal'
                        s_regen = phys.models[source_name].regenerate()
                        phys.models[source_name]['x'] = x
                        phys.models[source_name]['return_rate'] = return_rate
                        phys.models[source_name]['regen_mode'] = regen_mode
                        map_pores = phys.map_pores()
                        loc = pores[sp.in1d(pores, map_pores)]
                        if mode == 'merge':
                            try:
                                spore = self.pores('source_' + prop)
                                if sp.sum(sp.in1d(loc, spore)) > 0:
                                    raise Exception('Because of the existing '
                                                    'source term, the method '
                                                    'cannot apply new source '
                                                    'terms with the merge mode'
                                                    ' to the specified pores.')
                            except KeyError:
                                pass
                        self['pore.source_' + prop][loc] = True
                        map_pores_loc = sp.in1d(map_pores, pores)
                        self['pore.source_' + source_mode +
                             '_s1_' + prop][loc] = s_regen[:, 0][map_pores_loc]
                        self['pore.source_' + source_mode +
                             '_s2_' + prop][loc] = s_regen[:, 1][map_pores_loc]
                        if source_mode is not 'linear':
                            self['pore.source_maxiter'][loc] = maxiter
                            self['pore.source_tol'][loc] = tol
        else:
            Exception('No source_name has been sent for set_source_term ' +
                      'method in the algorithm ' + self.name)

    def run(self, **kwargs):
        r"""
        This calls the setup method in the algorithm and then runs the outer
        iteration stage.
        All of the arguments used in setup and solve methods, can be sent here
        as kwargs.
        """
        logger.info("Setup " + self.__class__.__name__)
        self.setup(**kwargs)
        self._do_outer_iteration_stage(**kwargs)

    def _do_outer_iteration_stage(self, **kwargs):
        r"""
        This calls the solve method in the algorithm.
        Many other outer loops can be added here as well, before or after
        calling solve method.
        """
        self.solve(**kwargs)

    def solve(self, A=None, b=None, iterative_solver=None, **kwargs):
        r"""
        Executes the right algorithm for the solution: regular solution of a
        linear system or iterative solution over the nonlinear source terms.

        Parameters
        ----------
        A : sparse matrix
            2D Coefficient matrix
        b : dense matrix
            1D RHS vector
        iterative_sovler : string
            Name of solver to use.  If not solve is specified, sp.solve is used
            which is a direct solver (SuperLU on default Scipy installation)
        kwargs : list of keyword arguments
            These arguments and values are sent to the sparse solver, so read
            the specific documentation for the solver chosen
        """
        self._iterative_solver = iterative_solver

        # Executes the right algorithm
        if any('pore.source_nonlinear' in s for s in self.props()):
            X = self._do_one_outer_iteration(**kwargs)
        else:
            X = self._do_one_inner_iteration(A, b, **kwargs)
        self.X = X
        self._Neumann_super_X = self.X[self.Np:self._coeff_dimension]
        # Removing the additional super pore variables from the results
        self[self._quantity] = self.X[self.Ps]
        logger.info('Writing the results to ' + '[\'' + self._quantity +
                    '\'] in the ' + self.name + ' algorithm.')

    def _do_one_inner_iteration(self, A, b, **kwargs):
        r"""
        This method solves AX = b and returns the result to the corresponding
        algorithm.
        """
        logger.info('Solving AX = b for the sparse matrices')

        if A is None:
            A = self.A
        if b is None:
            b = self.b
        if self._iterative_solver is None:
            X = sprslin.spsolve(A, b)
        else:
            params = kwargs.copy()
            solver_params = ['x0', 'tol', 'maxiter', 'xtype', 'M', 'callback']
            [params.pop(item, None) for item in kwargs.keys()
             if item not in solver_params]
            tol = kwargs.get('tol')
            if tol is None:
                tol = 1e-20
            params['tol'] = tol
            if self._iterative_solver == 'cg':
                result = sprslin.cg(A, b, **params)
            elif self._iterative_solver == 'gmres':
                result = sprslin.gmres(A, b, **params)
            elif self._iterative_solver == 'bicgstab':
                result = sprslin.bicgstab(A, b, **params)
            X = result[0]
            self._iterative_solver_info = result[1]
        return X

    def _do_one_outer_iteration(self, **kwargs):
        r"""
        One iteration of an outer iteration loop for an algorithm
        (e.g. time or parametric study)
        """
        # Checking for the necessary values in Picard algorithm
        nan_tol = sp.isnan(self['pore.source_tol'])
        nan_max = sp.isnan(self['pore.source_maxiter'])
        self._tol_for_all = sp.amin(self['pore.source_tol'][~nan_tol])
        self._maxiter_for_all = sp.amax(self['pore.source_maxiter'][~nan_max])
        if self._guess is None:
            self._guess = sp.zeros(self._coeff_dimension)
        t = 1
        step = 0
        # The main Picard loop
        while t > self._tol_for_all and step <= self._maxiter_for_all:
            X, t, A, b = self._do_inner_iteration_stage(guess=self._guess,
                                                        **kwargs)
            logger.info('tol for Picard source_algorithm in step ' +
                        str(step) + ' : ' + str(t))
            self._guess = X
            step += 1
        # Check for divergence
        self._steps = step
        if t >= self._tol_for_all and step > self._maxiter_for_all:
            raise Exception('Iterative algorithm for the source term reached '
                            'to the maxiter: ' + str(self._maxiter_for_all) +
                            ' without achieving tol: ' +
                            str(self._tol_for_all))
        logger.info('Picard algorithm for source term converged!')
        self.A = A
        self.b = b
        self._tol_reached = t
        return X

    def _do_inner_iteration_stage(self, guess, **kwargs):
        r"""
        This inner loop updates the source terms based on the new values of
        the quantity, then modifies A and b matrices, solves AX = b and
        returns the result.
        """
        # Updating the source terms
        s1 = sp.zeros(self._coeff_dimension)
        s2 = sp.zeros(self._coeff_dimension)
        for label in self.labels():
            if 'pore.source_' in label:
                source_name = label.replace('pore.source_', '')
                if 'pore.source_nonlinear_s1_' + source_name in self.props():
                    arr = self.pores('source_'+source_name)
                    tol = min(sp.unique(self['pore.source_tol'][arr]))
                    maxiter = max(sp.unique(self['pore.source_maxiter'][arr]))
                    self.set_source_term(source_name=source_name,
                                         pores=self.pores(label), x0=guess,
                                         tol=tol, maxiter=maxiter,
                                         mode='update')
                    prop1 = 'pore.source_nonlinear_s1_' + source_name
                    mask1 = ~sp.isnan(self[prop1])
                    s1[~sp.isnan(self[prop1])] = s1[mask1] + self[prop1][mask1]
                    prop2 = 'pore.source_nonlinear_s2_' + source_name
                    mask2 = ~sp.isnan(self[prop2])
                    s2[~sp.isnan(self[prop2])] = s2[mask2] + self[prop2][mask2]
        self.s1 = s1
        self.s2 = s2
        # Modifying A and b
        pores = self.pores('source_*')
        S1 = s1[pores]
        S2 = s2[pores]
        A = self._build_coefficient_matrix(modified_diag_pores=pores,
                                           diag_added_data=S1,
                                           mode='modify_diagonal')
        b = self._build_RHS_matrix(modified_RHS_pores=pores,
                                   RHS_added_data=-S2, mode='modify_RHS')
        # Solving AX = b
        X = self._do_one_inner_iteration(A=A, b=b, **kwargs)
        # Calculates absolute error
        t = sp.amax(sp.absolute(guess - X))
        return X, t, A, b

    def return_results(self, pores=None, throats=None, **kwargs):
        r"""
        Send results of simulation out the the appropriate locations.

        This is a basic version of the update that simply sends out the main
        result (quantity). More elaborate updates should be subclassed.
        """
        if pores is None:
            pores = self.Ps
        if throats is None:
            throats = self.Ts

        phase_quantity = self._quantity.replace(self._phase.name + '_', '')
        if phase_quantity not in self._phase.props():
            self._phase[phase_quantity] = sp.nan
        self._phase[phase_quantity][pores] = self[self._quantity][pores]
        conn_arr = self._net.find_connected_pores(self.Ts)
        dx = sp.squeeze(sp.diff(self[self._quantity][conn_arr], n=1, axis=1))
        g = self['throat.conductance']
        rate = sp.absolute(g * dx)
        if 'throat.rate' not in self._phase.props():
            self._phase['throat.rate'] = sp.nan
        self._phase['throat.rate'][throats] = rate[throats]
        logger.debug('Results of ' + self.name +
                     ' algorithm have been added to ' + self._phase.name)

    def _build_coefficient_matrix(self, modified_diag_pores=None,
                                  diag_added_data=None, mode='overwrite'):
        r"""
        This builds the sparse coefficient matrix for the linear solver.
        """
        if mode == 'overwrite':
            # Filling coefficient matrix
            tpore1 = self._net['throat.conns'][:, 0]
            tpore2 = self._net['throat.conns'][:, 1]
            # Identify Dirichlet pores
            try:
                temp = self.pores(self._phase.name + '_Dirichlet',
                                  mode='difference')
            except KeyError:
                temp = self.Ps
                logger.warning('No direct Dirichlet boundary condition has ' +
                               'been applied to the phase ' +
                               self._phase.name + ' in the algorithm ' +
                               self.name)
            loc1 = sp.in1d(tpore1, temp)
            loc2 = sp.in1d(tpore2, temp)
            modified_tpore1 = tpore1[loc1]
            modified_tpore2 = tpore2[loc1]
            row = modified_tpore1
            col = modified_tpore2
            # Expand the conductance to a vector if necessary
            g = self['throat.conductance']
            if sp.size(g) == 1:
                g = g * sp.ones(self.Nt)
            data_main = g
            data = data_main[loc1]
            modified_tpore2 = tpore2[loc2]
            modified_tpore1 = tpore1[loc2]
            row = sp.append(row, modified_tpore2)
            col = sp.append(col, modified_tpore1)
            data = sp.append(data, data_main[loc2])
            A_dim = self.Np
            # Check for Neuman_group BCs and add superpores if necessary
            try:
                self.pores(self._phase.name + '_Neumann_group')
                self._extra_Neumann_size = len(getattr(self, '_pore_' +
                                                       self._phase.name +
                                                       '_Neumann_group_' +
                                                       'location'))
                self._group_Neumann_vals = sp.zeros(self._extra_Neumann_size)
                l_g_super = len(self.super_pore_conductance)
                if l_g_super not in [0, 1, self._extra_Neumann_size]:
                    raise Exception('length of the list of super_pore_'
                                    'conductance and the number of different'
                                    ' Neumann_group BCs do not match.')
                if l_g_super == 1:
                    t = [sp.array(self.super_pore_conductance)]
                    self.super_pore_conductance = t * self._extra_Neumann_size
                for N in sp.arange(0, self._extra_Neumann_size):
                    neu_tpore2 = getattr(self, '_pore_' + self._phase.name +
                                         '_Neumann_group_location')[N]
                    Nval = self['pore.' + self._phase.name +
                                '_bcval_Neumann_group']
                    self._group_Neumann_vals[N] = sp.unique(Nval[neu_tpore2])
                    nt = self._net.find_neighbor_throats(pores=neu_tpore2)
                    try:
                        g_super = self.super_pore_conductance[N]
                    except IndexError:
                        g_super = 1e-3 * min(data_main[nt])
                        self.super_pore_conductance.append(g_super)
                    if sp.size(g_super) == 1:
                        g_super = len(neu_tpore2)*[g_super]
                    row = sp.append(row, neu_tpore2)
                    col = sp.append(col, len(neu_tpore2) * [A_dim + N])
                    data = sp.append(data, g_super)
                    row = sp.append(row, len(neu_tpore2) * [A_dim + N])
                    col = sp.append(col, neu_tpore2)
                    data = sp.append(data, g_super)
                A_dim = A_dim + self._extra_Neumann_size
            except KeyError:
                pass
            # Adding positions for diagonal
            diag = sp.arange(0, A_dim)
            try:
                pores = self.pores(self._phase.name + '_Dirichlet')
                row = sp.append(row, diag[pores])
                col = sp.append(col, diag[pores])
                data = sp.append(data, sp.ones_like(diag[pores]))
                temp_data = sp.copy(data)
                temp_data[sp.in1d(row, diag[pores])] = 0
                non_Dir_diag = diag[~sp.in1d(diag, diag[pores])]
            except KeyError:
                temp_data = sp.copy(data)
                non_Dir_diag = diag
            S_temp = sp.zeros(A_dim)
            for i in sp.arange(0, len(row)):
                S_temp[row[i]] = S_temp[row[i]] - temp_data[i]
            # Store values for modifying the diagonal in mode='modify_diagonal'
            self._non_source_row = row
            self._non_source_col = col
            self._non_source_data = data
            self._non_Dir_diag = non_Dir_diag
            self._diagonal_vals = S_temp
            self._coeff_dimension = A_dim

        if mode in ['overwrite', 'modify_diagonal']:
            diagonal_vals = sp.copy(self._diagonal_vals)
            # Adding necessary terms to the diagonal such as source terms
            if modified_diag_pores is not None and diag_added_data is not None:
                if sp.size(modified_diag_pores) == sp.size(diag_added_data):
                    sec1 = self._diagonal_vals[modified_diag_pores]
                    sec2 = diag_added_data
                    diagonal_vals[modified_diag_pores] = sec1 + sec2
                else:
                    raise Exception('Provided data and pores for modifying '
                                    'coefficient matrix should have the same' +
                                    ' size!')
                if mode == 'overwrite':
                    self._diagonal_vals = diagonal_vals
            data = sp.append(self._non_source_data,
                             diagonal_vals[self._non_Dir_diag])
            row = sp.append(self._non_source_row, self._non_Dir_diag)
            col = sp.append(self._non_source_col, self._non_Dir_diag)
            # Convert the lists to the sparse matrix
            a = sprs.coo.coo_matrix((data, (row, col)),
                                    (self._coeff_dimension,
                                     self._coeff_dimension))
            A = a.tocsr()
            A.eliminate_zeros()
            return(A)

    def _build_RHS_matrix(self, modified_RHS_pores=None, RHS_added_data=None,
                          mode='overwrite'):
        r"""
        This builds the right-hand-side matrix for the linear solver.
        """

        if mode == 'overwrite':
            A_dim = self._coeff_dimension
            b = sp.zeros([A_dim, 1])
            try:
                Dir_pores = self.pores(self._phase.name + '_Dirichlet')
                Dir_pores_vals = self['pore.' + self._phase.name +
                                      '_bcval_Dirichlet'][Dir_pores]
                b[Dir_pores] = sp.reshape(Dir_pores_vals, [len(Dir_pores), 1])
            except KeyError:
                pass
            try:
                ind_Neu_pores = self.pores(self._phase.name + '_Neumann')
                ind_Neu_pores_vals = self['pore.' + self._phase.name +
                                          '_bcval_' + 'Neumann'][ind_Neu_pores]
                b[ind_Neu_pores] = sp.reshape(ind_Neu_pores_vals,
                                              [len(ind_Neu_pores), 1])
            except KeyError:
                pass
            try:
                self.pores(self._phase.name + '_Neumann_group')
                pnum = self._net.Np
                NG_loc = sp.r_[pnum: (pnum + len(self._group_Neumann_vals))]
                NG_l = len(self._group_Neumann_vals)
                NG_arr = self._group_Neumann_vals[sp.r_[0:NG_l]]
                b[NG_loc] = sp.reshape(NG_arr, [NG_l, 1])
            except KeyError:
                pass
        if mode in ['overwrite', 'modify_RHS']:
            try:
                b = sp.copy(self.b)
            except AttributeError:
                pass
            # Adding necessary terms to the RHS for non-Dirichlet pores
            if modified_RHS_pores is not None and RHS_added_data is not None:
                if sp.size(modified_RHS_pores) == sp.size(RHS_added_data):
                    p = sp.in1d(modified_RHS_pores, self._non_Dir_diag)
                    data = RHS_added_data[p]
                    b[modified_RHS_pores[p]] = b[modified_RHS_pores[p]] + \
                                               data.reshape([len(data), 1])
                else:
                    raise Exception('Provided data and pores for modifying'
                                    ' RHS matrix should have the same size!')
        return(b)

    def rate(self, pores=None, network=None, conductance=None, X_value=None,
             mode='group'):
        r"""
        Send a list of pores and receive the net rate
        of material moving into them.

        Parameters
        ----------
        pores : array_like
            The pores where the net rate will be calculated
        network : OpenPNM Network Object
            The network object to which this algorithm will apply.
            If no network is sent, the rate will apply to the network which is
            attached to the algorithm.
        conductance : array_like
            The conductance which this algorithm will use to calculate the
            rate.
            If no conductance is sent, the rate will use the
            'throat.conductance' which is attached to the algorithm.
        X_value : array_like
            The values of the quantity (temperature, mole_fraction,
            voltage, ...), which this algorithm will use to calculate the rate.
            If no X_value is sent, the rate will look at the '_quantity',
            which is attached to the algorithm.
        mode : string, optional
            Controls how to return the rate.  Options are:
            - 'group'(default): It returns the cumulative rate moving into them
            - 'single': It calculates the rate for each pore individually.
        """

        if network is None:
            network = self._net
        if conductance is None:
            conductance = self['throat.conductance']
        if X_value is None:
            X_value = self[self._quantity]
        pores = sp.array(pores, ndmin=1)
        R = []
        if mode == 'group':
            t = network.find_neighbor_throats(pores, flatten=True,
                                              mode='not_intersection')
            throat_group_num = 1
        elif mode == 'single':
            t = network.find_neighbor_throats(pores, flatten=False,
                                              mode='not_intersection')
            throat_group_num = sp.size(t)
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
            X1 = X_value[pores1]
            X2 = X_value[pores2]
            g = conductance[throats]
            R.append(sp.sum(sp.multiply(g, (X2 - X1))))
        return(sp.array(R, ndmin=1))

    def _calc_eff_prop(self, check_health=False):
        r"""
        This returns the main parameters for calculating the effective
        property in a linear transport equation.
        It also checks for the proper boundary conditions, inlets and outlets.

        Parameters
        ----------
        check_health : boolean(optional)
            It analyzes the inlet and outlet pores to check their spatial
            positions
        """
        try:
            self[self._quantity]
        except KeyError:
            raise Exception('The algorithm has not been run yet. Cannot ' +
                            'calculate effective property.')
        # Determine boundary conditions by analyzing algorithm object
        Ps = self.pores('pore.' + self._phase.name + '_Dirichlet')
        BCs = sp.unique(self['pore.' + self._phase.name +
                             '_bcval_Dirichlet'][Ps])
        if sp.shape(BCs)[0] != 2:
            raise Exception('The supplied algorithm did not have appropriate' +
                            ' BCs')
        inlets = sp.where(self['pore.' + self._phase.name +
                               '_bcval_Dirichlet'] == sp.amax(BCs))[0]
        outlets = sp.where(self['pore.' + self._phase.name +
                                '_bcval_Dirichlet'] == sp.amin(BCs))[0]

        # Analyze input and output pores
        if check_health:
            # Check for coplanarity
            if self._net.iscoplanar(inlets) == False:
                raise Exception('The inlet pores do not define a plane. ' +
                                'Effective property will be approximation')
            if self._net.iscoplanar(outlets) == False:
                raise Exception('The outlet pores do not define a plane. ' +
                                'Effective property will be approximation')
            # Ensure pores are on a face of domain
            # (only 1 non-self neighbor each)
            PnI = self._net.find_neighbor_pores(pores=inlets,
                                                mode='not_intersection',
                                                excl_self=True)
            if sp.shape(PnI) != sp.shape(inlets):
                logger.warning('The inlet pores have too many neighbors. ' +
                               'Internal pores appear to be selected.')
                pass
            PnO = self._net.find_neighbor_pores(pores=outlets,
                                                mode='not_intersection',
                                                excl_self=True)
            if sp.shape(PnO) != sp.shape(outlets):
                logger.warning('The outlet pores have too many neighbors. ' +
                               'Internal pores appear to be selected.')
                pass

        # Fetch area and length of domain
        if 'pore.vert_index' in self._net.props():
            A = vo.vertex_dimension(network=self._net, face1=inlets,
                                    parm='area')
            L = vo.vertex_dimension(network=self._net, face1=inlets,
                                    face2=outlets, parm='length')
        else:
            A = self._net.domain_area(face=inlets)
            L = self._net.domain_length(face_1=inlets, face_2=outlets)
        flow = self.rate(pores=inlets)
        D = sp.sum(flow)*L/A/(BCs[0] - BCs[1])
        return D
