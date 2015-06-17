# -*- coding: utf-8 -*-
"""
===============================================================================
module __StokesFlow__: Viscous fluid flow
===============================================================================

"""
import scipy as sp
from OpenPNM.Algorithms import GenericLinearTransport
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class StokesFlow(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the network.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,
    ...                                     pores=pn.pores(),
    ...                                     throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,
    ...                                     pores=pn.pores(),throats=pn.throats())
    >>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=phase1)
    >>> BC1_pores = pn.pores('top')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    >>> BC2_pores = pn.pores('bottom')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    >>> alg.run()
    >>> alg.return_results()
    >>> Peff = round(alg.calc_eff_permeability(), 10)
    >>> print(Peff)
    1.8663e-05
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')

    def setup(self, conductance='hydraulic_conductance', quantity='pressure',
              super_pore_conductance=None, **params):
        r"""
        This setup provides the initial requirements for the solver setup.
        """
        logger.info('Setup ' + self.__class__.__name__)
        super().setup(conductance=conductance, quantity=quantity,
                      super_pore_conductance=super_pore_conductance)

    def calc_eff_permeability(self):
        r"""
        This calculates the effective permeability in this linear
        transport algorithm.
        """
        d_normal = self._calc_eff_prop()
        self._eff_property = d_normal * sp.mean(self._phase['pore.viscosity'])
        return self._eff_property
