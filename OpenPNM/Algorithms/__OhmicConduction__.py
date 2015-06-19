# -*- coding: utf-8 -*-
"""
===============================================================================
module __OhmicConduction__: Electronic or ionic conduction
===============================================================================

"""
import scipy as sp
from OpenPNM.Algorithms import GenericLinearTransport
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class OhmicConduction(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate electron and ionic
    conduction.  The 2 main roles of this subclass are to set the default
    property names and to implement a method for calculating the effective
    conductivity of the network.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn, pores=pn.pores(),
    ...                                     throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,
    ...                                     pores=pn.pores(), throats=pn.throats())
    >>> phys1['throat.electrical_conductance'] = 1
    >>> alg = OpenPNM.Algorithms.OhmicConduction(network=pn, phase=phase1)
    >>> BC1_pores = pn.pores('top')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    >>> BC2_pores = pn.pores('bottom')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    >>> alg.run()
    >>> alg.return_results()
    >>> Ceff = round(alg.calc_effective_conductivity(), 3)
    >>> print(Ceff)
    1.012

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')

    def setup(self, conductance='electrical_conductance',
              quantity='voltage', super_pore_conductance=None, **params):
        r"""
        This setup provides the initial requirements for the solver setup.
        """
        logger.info('Setup ' + self.__class__.__name__)
        super().setup(conductance=conductance, quantity=quantity,
                      super_pore_conductance=super_pore_conductance)

    def calc_effective_conductivity(self):
        r"""
        This calculates the effective electrical conductivity in this linear
        transport algorithm.
        """
        return self._calc_eff_prop()
