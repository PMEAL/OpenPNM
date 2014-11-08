# -*- coding: utf-8 -*-
"""
===============================================================================
module __FourierConduction__: Conductive heat transfer
===============================================================================

A subclass of GenericLinearTransport to simulate heat conduction

"""
import scipy as sp
from OpenPNM.Algorithms import GenericLinearTransport
from OpenPNM.Base import logging
logger = logging.getLogger()

class FourierConduction(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the effective conductivity of the network.

    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
    >>> alg = OpenPNM.Algorithms.FourierConduction(network=pn, phase=phase1)
    >>> BC1_pores = pn.pores('top')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    >>> BC2_pores = pn.pores('bottom')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    >>> alg.run()
    >>> alg.return_results()
    >>> Ceff = round(alg._calc_eff_prop(), 3) #This line and the next line should fail until someone writes this function
    >>> print(Ceff) #unless something changed with our test objects, this should print "0.025"
    0.822


    """

    def __init__(self,**kwargs):
        r'''
        '''
        super(FourierConduction,self).__init__(**kwargs)
        logger.info('Create '+self.__class__.__name__+' Object')

    def run(self,conductance='thermal_conductance',quantity='temperature',**params):
        r'''
        '''
        logger.info('Setup '+self.__class__.__name__)
        super(FourierConduction,self).setup(conductance=conductance,quantity=quantity)

        super(GenericLinearTransport,self).run()


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
