# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:50:28 2013

@author: Mahmoudreza Aghighi

module __FickianDiffusion__: Fick's Law Diffusion
========================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.sparse as sprs
import scipy.sparse.linalg as splin
from __TransportsLinearSolver__ import TransportsLinearSolver

class FickianDiffusion(TransportsLinearSolver):
    r"""   
    
    FickianDiffusion - Class to run Fick's law mass transfer diffusion on constructed networks
    
                        It returns conecentration gradient inside the network.
                        An invasion algorithm should be used before running diffusion calculations.
                                   
                            
    Parameters
    ----------
    - Alg:
        Algorithm for Non wetting phase configuration
        'OP' , 'IP' and 'None' are possible options.        
    - Pressure for 'OP':
        Applied pressure on the network which causes some pores and throats be invaded by non-wetting phase.
        It can be a single value or a list.
        The class will consider uninvaded conduits based on this pressure, as possible voids for gas diffusion.
        Every pore with Pc_invade greater than Pressure, will be considered as open pores        
    - Psequence for 'IP':
        The desired sequence of invaded pores based on IP algorithm.
        The class will consider uninvaded conduits based on this sequence number, as possible voids for gas diffusion.
        It can be a single value or a list.
        Every pore with sequence number greater than Psequence, will be considered as open pores.
    - Pinvaded for 'None':
        It is a list which contains invaded pores and based on user's definition.
            
    -loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)    
    - Total_Conc: 
        Total gas concentration
    - Diff_Coefficient:
        Diffision coefficient for diffusion of A through stagnant B(e.g. Oxygen through Nitrogen and Water vapor)
                
    Examples
    --------

    All of the variables in this class have default values, but users can define them too.
    >>>
    
    To Do:
        - Instead of sending (Pressure in OP) or (Psequence in IP) as inputs, this class should only accept 
        a list of invaded voids.
        - loglevel and documentation is not complete yet.
        - For Applying Neumann boundary condition, a function is needed to consider mass conservation, sink or
        source terms and warn user about the applied values.
        
    """
    
    def __init__(self,net,loglevel=10,**kwargs):
        r"""
        Initializing the class
        """
        super(FickianDiffusion,self).__init__(net = net,**kwargs)
        self._logger.info("Create Fick's Diffusion Algorithm Object")
        print 'init'
            
    def _setup_for_TransportSolver(self,                                    
                                    loglevel=10,
                                    DryNetwork=0,
                                    T=353.15,
                                    P=1.01325e5,                                    
                                    **params):
        r"""
        Main features: Applying Boundary Conditions & Creating Transport Conductivities 

        This function executes the essential mathods for building matrices in Linear solution 
        """
        print '_setup_for_FickianDiffusion'
        self._net.pore_properties['temperature'] = T
        self._net.pore_properties['GasPressure'] = P
        self.DryNetwork = DryNetwork
        if 'Cdiff' not in self._net.throat_properties:
            self._logger.info("Creating diffusion conductivities for dry conduits")
            self._net.throat_properties['Cdiff_dry'] = OpenPNM.Physics.MassTransport.Conduits_DiffusionConductivity(self._net,**params) 
        else:
            self.logger.info("Diffusion conductivities for dry conduits have already been created")     

        if not self.DryNetwork:
            self.logger.info("Applying non-wetting phase saturation states to conduit conductivities")
            OpenPNM.Physics.MultiPhase.Conduit_Filled_State_Calculator(self._net,**params)
            self._net.throat_properties['Conductivities_Exp'] = OpenPNM.Physics.MultiPhase.Apply_Phase_State_to_Conduit_Conductivity(self._net,**params)
        else:
            self._net.throat_properties['Conductivities_Exp'] = self._net.throat_properties['Cdiff_dry']

    
    def _do_inner_iteration_stage(self,Experiment='Exp1'):
        r"""
                       
        """
        print '_do_outer_iteration_stage'
        self._net.pore_properties[Experiment+'_mole_fractions'] = self._do_one_inner_iteration()

