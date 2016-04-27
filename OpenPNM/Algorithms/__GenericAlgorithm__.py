# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
==================================================================

This generic class contains the recommended methods for subclassed algorithms.
It inherits from Core, so is Python Dict with the OpenPNM data control methods.

"""
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
            network = GenericNetwork()
        self.network.update({network.name: network})

        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']

    def set_boundary_conditions(self, **kwargs):
        r"""
        Not implemented
        """
        raise NotImplementedError('This method must be implemeted by each ' +
                                  'specific Algorithm class')

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
        raise NotImplementedError('This method must be implemeted by each ' +
                                  'specific Algorithm class')

    def return_results(self, **kwargs):
        r"""
        Not implemented
        """
        raise NotImplementedError('This method must be implemeted by each ' +
                                  'specific Algorithm class')
