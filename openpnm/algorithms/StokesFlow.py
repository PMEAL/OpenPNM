# -*- coding: utf-8 -*-
"""
===============================================================================
module __StokesFlow__: Viscous fluid flow
===============================================================================

"""

import scipy as sp
from openpnm.algorithms import GenericLinearTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class StokesFlow(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the network.

    Examples
    --------
    >>> import openpnm as op
    >>> net = openpnm.network.Cubic(shape=[5, 5, 5])
    >>> geom = op.geometry.StickAndBall(network=net,
    ...                                 pores=net.Ps,
    ...                                 throats=net.Ts)
    >>> water = op.phases.Water(network=net)
    >>> phys_water = op.physics.GenericPhysics(network=net,
    ...                                        phase=water,
    ...                                        geometry=geom)
    >>> sf = openpnm.algorithms.StokesFlow(network=net, phase=water)
    >>> sf.set_BC(pores=net.pores('left'), bctype='dirichlet', bcvalue=1.0)
    >>> sf.set_BC(pores=net.pores('reft'), bctype='dirichlet', bcvalue=0.0)
    >>> sf.run()
    >>> Ke = sf.calc_eff_permeability() # Effective permeability
    >>> print('Effective permeability:', Ke)
    1.8663e-05
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')

    def setup(self, conductance='hydraulic_conductance',
              quantity='pressure', **kwrgs):
        r"""
        This setup provides the initial requirements for the solver setup.
        """
        logger.info('Setup ' + self.__class__.__name__)
        super().setup(conductance=conductance, quantity=quantity)

    def calc_eff_permeability(self):
        r"""
        This calculates the effective permeability in this linear
        transport algorithm.
        """
        phase = self.project.phases[self['phase']]
        d_normal = self._calc_eff_prop()
        self._eff_property = d_normal / sp.mean(phase['pore.viscosity'])
        return self._eff_property
