
"""
module Diffusivity
===============================================================================

"""
import OpenPNM
import scipy as sp

def set_as(fluid=None,diff=2.09e-5):
    r"""
    """
    diff = sp.array(diff,ndmin=1)
    fluid.update({'diffusivity':diff})

def ChapmanEnskog(fluid,T,P):
    r"""
    nothing yet
    """
    print 'nothing yet'
    return

def Fuller(fluid1,T,P,MA=31.99,MB=28.01,VA=16.6,VB=17.9):
    r"""
    Uses the Fuller model estimate diffusion coeffient from first principles at conditions of interest

    Parameters
    ----------
    T & P :  float, array_like
        Temperature and pressure of interest
    """
    MAB = 2*(1/MA+1/MB)**(-1)
    DAB = 0.00143*T**1.75/(P*MAB**0.5*(VA**(1/3)+VB**(1/3))**2)*1e-4
    OpenPNM.Diffusivity.set_as(fluid1,DAB)

def FullerScaling(fluid,DABo=2.09e-5,To=298.,Po=101325.,Ti=298.,Pi=101325.):
    r"""
    Uses the Fuller model to adjust a diffusion coeffciient from reference conditions to conditions of interest

    Parameters
    ----------
    network : OpenPNM Network Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To, Pi & Ti : float, array_like
        Pressure & temperature at reference conditions, pressure & temperature at conditions of interest, respectively
    """
    Di = DABo*(Po/To**1.75)*(Ti**1.75/Pi)
    OpenPNM.Fluids.Diffusivity.set_as(fluid,Di)
