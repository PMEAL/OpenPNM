# -*- coding: utf-8 -*-
"""
===============================================================================
module __FickianDiffusion__: Diffusive mass transfer
===============================================================================

"""

import scipy as sp
from openpnm.algorithms import GenericLinearTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class FickianDiffusion(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate binary diffusion. The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the effective diffusion coefficient
    of the network.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def setup(self, conductance='diffusive_conductance',
              quantity='mole_fraction', **kwargs):
        r"""

        """
        super().setup(conductance=conductance, quantity=quantity)

    def calc_eff_diffusivity(self):
        r"""
        This calculates the effective diffusivity in this linear transport
        algorithm.
        """
        phase = self.project.phases[self['phase']]
        d_normal = self._calc_eff_prop()
        self._eff_property = d_normal / sp.mean(phase['pore.molar_density'])
        return self._eff_property
