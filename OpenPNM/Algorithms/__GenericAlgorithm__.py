# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
==================================================================

This generic class contains the recommended methods for subclassed algorithms.
It inherits from Core, so is Python Dict with the OpenPNM data control methods.

"""
import scipy as sp
from OpenPNM.Base import Core
from OpenPNM.Base import logging
from OpenPNM.Network import GenericNetwork
logger = logging.getLogger(__name__)


class GenericAlgorithm(Core):
    r"""
    GenericAlgorithm - Base class to execute algorithms

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object to which this algorithm will apply.

    name : string, optional
        Name of this algorithm


    Notes
    -----
    If no network is supplied an empty algorithm object is returned.  This is
    useful for loading in a saved algorithm from memory.

    """

    def __init__(self, network=None, **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network

        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']

    def run(self, **params):
        r"""
        Main run command for the algorithm
        """
        self._do_outer_iteration_stage()

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._do_one_outer_iteration()

    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm
        (e.g. time or parametric study)
        """
        self._do_inner_iteration_stage()

    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._do_one_inner_iteration()

    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        pass

    def return_results(self, **kwargs):
        pass

    def set_boundary_conditions(self, component=None, bctype='', bcvalue=None,
                                pores=None, throats=None, mode='merge'):
        r"""
        Apply boundary conditions to specified pores or throats

        Parameters
        ----------
        bctype : string
            Specifies the type or the name of boundary condition to apply. \
            The types can be one one of the followings:
                 - 'Dirichlet' : Specify the quantity in each location
                 - 'Neumann' : Specify the flow rate into each location
                 - 'Neumann_group' : Specify the net flow rate into a group
                   of pores/throats
        component : OpenPNM Phase object
            The Phase object to which this BC applies
        bcvalue : array_like
            The boundary value to apply, such as concentration or rate
        pores : array_like
            The pores where the boundary conditions should be applied
        throats : array_like
            The throats where the boundary conditions should be applied
        mode : string, optional
            Controls how the conditions are applied.  Options are:

            - 'merge': Inserts the specified values, leaving existing values \
              elsewhere
            - 'overwrite': Inserts specified values, clearing all other values
            - 'remove': Removes boundary conditions from specified locations

        Notes
        -----
        - It is not possible to have multiple boundary conditions for a
          specified location in just one algorithm. So when new condition is
          going to be applied to a specific location, any existing one should
          be removed or overwritten.
        - BCs for pores and for throats should be applied independently.
        """
        try:
            self._existing_BC
        except AttributeError:
            self._existing_BC = []
        if component is None:
            if sp.size(self._phases) != 1:
                raise Exception('In each use of set_boundary_conditions ' +
                                'method, one component should be specified ' +
                                'or attached to the algorithm.')
            else:
                component = self._phases[0]
        else:
            if sp.size(component) != 1:
                raise Exception('For using set_boundary_conditions method, ' +
                                'only one component should be specified.')

        if mode not in ['merge', 'overwrite', 'remove']:
            raise Exception('The mode (' + mode + ') cannot be applied to ' +
                            'the set_boundary_conditions!')

        logger.debug('BC applies to the component: ' + component.name)
        # If mode is 'remove', also bypass checks
        if mode == 'remove':
            if pores is None and throats is None:
                if bctype == '':
                    raise Exception('No bctype/pore/throat is specified')
                else:
                    for item in self.labels():
                        item_spl = item.split('.')
                        if bctype == (item_spl[-1]).replace(self._phase.name +
                                                            '_', ''):
                            element = item_spl[0]
                            try:
                                del self[element + '.' + component.name +
                                         '_bcval_' + bctype]
                            except KeyError:
                                pass
                            try:
                                del self[element + '.' + component.name +
                                         '_' + bctype]
                            except KeyError:
                                pass
                    logger.debug('Removing ' + bctype + ' from all locations' +
                                 ' for ' + component.name + ' in ' +
                                 self.name)
                    self._existing_BC.remove(bctype)
            else:
                if pores is not None:
                    if bctype != '':
                        prop_label = 'pore.' + component.name + '_bcval_'\
                                     + bctype
                        self[prop_label][pores] = sp.nan
                        info_label = 'pore.' + component.name + '_' + bctype
                        self[info_label][pores] = False
                        logger.debug('Removing ' + bctype + ' from the ' +
                                     'specified pores for ' + component.name +
                                     ' in ' + self.name)
                    else:
                        raise Exception('Cannot remove BC from the pores ' +
                                        'unless bctype is specified')

                if throats is not None:
                    if bctype != '':
                        prop_label = 'throat.' + component.name + '_bcval_'\
                                     + bctype
                        self[prop_label][throats] = sp.nan
                        info_label = 'throat.' + component.name + '_' + bctype
                        self[info_label][throats] = False
                        logger.debug('Removing ' + bctype + ' from the ' +
                                     'specified throats for ' +
                                     component.name + ' in ' + self.name)
                    else:
                        raise Exception('Cannot remove BC from the throats ' +
                                        'unless bctype is specified')

            return
        # Validate bctype
        if bctype == '':
            raise Exception('bctype must be specified')
        # Validate pores/throats
        if pores is None and throats is None:
            raise Exception('pores/throats must be specified')
        elif pores is not None and throats is not None:
            raise Exception('BC for pores and throats must be specified ' +
                            'independently.')
        elif throats is None:
            element = 'pore'
            loc = sp.array(pores, ndmin=1)
            all_length = self.Np
        elif pores is None:
            element = 'throat'
            loc = sp.array(throats, ndmin=1)
            all_length = self.Nt
        else:
            raise Exception('Problem with the pore and/or throat list')
        # Validate bcvalue
        if bcvalue is not None:
            # Check bcvalues are compatible with bctypes
            if bctype == 'Neumann_group':  # Only scalars are acceptable
                if sp.size(bcvalue) != 1:
                    raise Exception('When specifying Neumann_group, bcval ' +
                                    'should be a scalar')
                else:
                    bcvalue = sp.float64(bcvalue)
                    if 'Neumann_group' not in self._existing_BC:
                        setattr(self, '_' + element + '_' + component.name +
                                '_Neumann_group_location', [])
                    getattr(self, '_' + element + '_' + component.name +
                            '_Neumann_group_location').append(loc)
            else:  # Only scalars or Np/Nt-long are acceptable
                if sp.size(bcvalue) == 1:
                    bcvalue = sp.ones(sp.shape(loc)) * bcvalue
                elif sp.size(bcvalue) != sp.size(loc):
                    raise Exception('The pore/throat list and bcvalue list ' +
                                    'are different lengths')
        # Confirm that prop and label arrays exist
        l_prop = element + '.' + component.name + '_bcval_' + bctype
        if l_prop not in self.props():
            self[l_prop] = sp.ones((all_length,), dtype=float) * sp.nan
        l_label = element + '.' + component.name + '_' + bctype
        if l_label not in self.labels():
            self[l_label] = sp.zeros((all_length,), dtype=bool)
        # Check all BC from specified locations, prior to setting new ones
        for item in self.labels():
            bcname = (item.split('.')[-1]).replace(component.name + '_', "")
            if bcname in self._existing_BC and item.split('.')[0] == element:
                if mode == 'merge':
                    try:
                        c1 = element + '.' + component.name
                        c2 = '_bcval_' + bcname
                        c1_label = c1 + c2
                        self[c1_label][loc]
                        condition1 = sp.isnan(self[c1_label][loc]).all()
                        c2_label = c1 + '_' + bcname
                        condition2 = sp.sum(self[c2_label][loc]) == 0
                        if not (condition1 and condition2):
                            raise Exception('Because of the existing BCs, ' +
                                            'the method cannot apply new BC ' +
                                            'with the merge mode to the ' +
                                            'specified pore/throat.')
                    except KeyError:
                        pass
        # Set boundary conditions based on supplied mode
        if mode == 'merge':
            if bcvalue is not None:
                self[l_prop][loc] = bcvalue
            self[l_label][loc] = True
            if bctype not in self._existing_BC:
                self._existing_BC.append(bctype)
        elif mode == 'overwrite':
            self[l_prop] = sp.ones((all_length,), dtype=float) * sp.nan
            if bcvalue is not None:
                self[l_prop][loc] = bcvalue
            self[l_label] = sp.zeros((all_length,), dtype=bool)
            self[l_label][loc] = True
            if bctype not in self._existing_BC:
                self._existing_BC.append(bctype)
