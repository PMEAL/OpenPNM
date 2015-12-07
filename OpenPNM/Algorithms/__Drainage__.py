# -*- coding: utf-8 -*-
"""
===============================================================================
module __Drainage__: Capillary Pressure Curve Simulation
===============================================================================

"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class PrimaryDrainage(GenericAlgorithm):
    r"""
    Simulates a capillary drainage experiment by applying a list of increasing
    capillary pressures and invading all throat that are accessible and
    invadable at the given pressure

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    name : string, optional
        The name to assign to the Algorithm Object

    Notes
    -----
    1. Trapping of the wetting phase can be included

    """

    def __init__(self, network, name=None, **kwargs):
        super().__init__(network=network, name=name)
        if len(kwargs.keys()) > 0:
            self.setup(**kwargs)

    def setup(self):
        r"""
        This method is used to specify necessary arguments to the simulation.

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.
        throat_prop : string
            The dictionary key on the Phase object where the throat entry
            pressure values can be found.
        trapping : boolean
            Specifies whether wetting phase trapping should be included or not.
            Note that wetting phase outlets can be specified using the
            ``set_outlets`` method.  Otherwise it is assumed the wetting phase
            has no outlets.

        Notes
        -----
        In older versions many of these arguments were passed directly to the
        initialization of the Algorithm.  This is still supported by
        automatically calling ``setup`` from ``init`` if extra arguments are
        received.
        """
        pass

    def set_inlets(self, pores):
        r"""
        """

    def set_outlets(self, pores):
        r"""
        """
        pass

    def run(self, npts=25, inv_points=None, access_limited=True, **kwargs):
        r"""
        """
        pass

    def _do_outer_iteration_stage(self, inv_points):
        r"""
        """
        pass

    def _do_one_inner_iteration(self, inv_val):
        r"""
        """
        pass