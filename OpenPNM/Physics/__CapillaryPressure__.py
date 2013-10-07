#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericPhysics__: Base class to define pore scale physics
==================================================================

.. warning:: The classes of this module should be loaded through the 'Physics/__init__.py' file.

"""

import OpenPNM
import scipy as sp

class CapillaryPressure(OpenPNM.Utilities.OpenPNMbase):
    r"""
    Methods in this class are used to determine the capillary entry pressure of throats from their geometric and material properties.  
    """
    
    def __init__(self,**kwargs):
        super(CapillaryPressure,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        
    def Washburn(self,network,sigma,theta):
        r"""
        Computes the capillary entry pressure assuming the throat is a cylindrical tube.
        
        Parameters
        ----------
        network : OpenPNM Network Object
            The network to apply the calculation 
        
        sigma : float
            Surface tension of the invading/defending fluid pair.  Units must be consistent with the throat size values, but SI is encouraged.
    
        theta : float
            Contact angle formed by a droplet of the invading fluid and solid surface, measured through the defending fluid phase.  Angle must be given in degrees.
        
        Notes
        -----
        The Washburn equation is:
        
        .. math::
            P_c = -\frac{2\sigma(cos(\theta))}{r}

        This is the most basic approach to calcualing entry pressure and is suitable for highly non-wetting invading fluids in most materials.  
        
        """
        self._logger.info("Calculating throat entry capillary pressure using the Washburn equation")
        return -4*sigma*sp.cos(sp.radians(theta))/network.throat_properties['diameter']
        
    def Purcell(self,network,sigma,theta,r_toroid):
        r"""
        Computes the throat capillary entry pressure assuming the throat is a toroid.
        
        Parameters
        ----------
        network : OpenPNM Network Object
            The network to apply the calculation 
        
        sigma : float
            Surface tension of the invading/defending fluid pair.  Units must be consistent with the throat size values, but SI is encouraged.
    
        theta : float
            Contact angle formed by a droplet of the invading fluid and solid surface, measured through the defending fluid phase.  Angle must be given in degrees.
        
        r_toroid : float or array_like
            The radius of the solid
            
        Notes
        -----
        This approach accounts for the converging-diverging nature of many throat types.  Advancing the meniscus beyond the apex of the toroid requires an increase in capillary pressure beyond that for a cylindical tube of the same radius. The details of this equation are described by Mason and Morrow [1]_, and explored by Gostick [2]_ in the context of a pore network model.

        References
        ----------
        
        .. [1] G. Mason, N. R. Morrow, Effect of contact angle on capillary displacement curvatures in pore throats formed by spheres. J. Colloid Interface Sci. 168, 130 (1994).
        .. [2] J. Gostick, Random pore network modeling of fibrous PEMFC gas diffusion media using Voronoi and Delaunay tessellations. J. Electrochem. Soc. 160, F731 (2013).
        
        """
        #This seesm to work, but I wrote it quickly and lost track of the degree-radians conversions
        """TODO: 
        Triple check the accuracy of this equation
        """
        r = network.throat_properties['diameter']/2
        R = r_toroid
        alpha = theta - 180 + sp.arcsin(sp.sin(sp.radians(theta)/(1+r/R)))
        Pc = (-2*sigma/r)*(sp.cos(sp.radians(theta - alpha))/(1 + R/r*(1-sp.cos(sp.radians(alpha)))))
        return Pc    
        
    def Morrow(self,network,sigma,theta):
        r"""
        Computes the throat capillary pressure using simplified version of the Purcell toroid
        
        Parameters
        ----------
        network : OpenPNM Network Object
            The network to apply the calculation 
        
        sigma : float
            Surface tension of the invading/defending fluid pair.  Units must be consistent with the throat size values, but SI is encouraged.
    
        theta : float
            Contact angle formed by a droplet of the invading fluid and solid surface, measured through the defending fluid phase.  Angle must be given in degrees.
            
        Notes
        -----
        Mason and Morrow compared the Purcell toroid to experimental data on various sized monodisperse PTFE beads.  They found that data could be approximated decently by simply scaling the contact angle measured through the wetting phase by 2/3.  
        """
        return -4*sigma*sp.cos(sp.radians(2/3*(180-theta)))/network.throat_properties['diameter']