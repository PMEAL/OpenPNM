"""
module __EffectiveProperty__: Base class to estimate transport properties
===============================================================================

"""
import OpenPNM
import scipy as sp
import scipy.signal as spsg
import scipy.spatial as sptl

from .__GenericAlgorithm__ import GenericAlgorithm

class EffectiveProperty(GenericAlgorithm):
    r'''
    '''
    def __init__(self,**kwargs):
        r'''
        '''
        super(EffectiveProperty,self).__init__(**kwargs)
        self._logger.info("Construct Algorithm")
        
        
    def setup(self,algorithm,fluid,conductance=str,quantity=str):
        r'''
        '''
        self._alg = algorithm
        self._fluid = fluid
        self._conductance = 'throat.'+conductance.split('.')[-1]
        self._quantity = 'pore.'+quantity.split('.')[-1]
        
    def run(self):
        r'''
        '''
        #Determine boundary conditions by analyzing algorithm object
        Ps = self._alg.pores(labels='pore.Dirichlet')
        BCs = sp.unique(self._alg['pore.bcval_Dirichlet'][Ps])
        if sp.shape(BCs)[0] != 2:
            raise Exception('The supplied algorithm did not have appropriate BCs')
        inlets = sp.where(self._alg['pore.bcval_Dirichlet']==sp.amax(BCs))[0]
        outlets = sp.where(self._alg['pore.bcval_Dirichlet']==sp.amin(BCs))[0]

        #Analyze input and output pores
        #Check for coplanarity
        if self._net.iscoplanar(inlets) == False:
            raise Exception('The inlet pores do not define a plane')
        if self._net.iscoplanar(outlets) == False:
            raise Exception('The outlet pores do not define a plane')
        #Ensure pores are on a face of domain (only 1 non-self neighbor each)
        PnI = self._net.find_neighbor_pores(pores=inlets,mode='not_intersection',excl_self=True)
        if sp.shape(PnI) != sp.shape(inlets):
            raise Exception('The inlet pores have too many neighbors')
        PnO = self._net.find_neighbor_pores(pores=outlets,mode='not_intersection',excl_self=True)
        if sp.shape(PnO) != sp.shape(outlets):
            raise Exception('The outlet pores have too many neighbors')
        Pin = inlets
        Pout = outlets
        
        #Fetch area and length of domain
        A = self._net.domain_area(face=Pin)
        L = self._net.domain_length(face_1=Pin,face_2=Pout)
        
        #Find flow through inlet face
        Pn = self._net.find_neighbor_pores(pores=Pin,excl_self=True)
        Ts = self._net.find_connecting_throat(Pin,Pn)
        g = self._fluid[self._conductance][Ts]
        s = self._fluid['throat.occupancy'][Ts]
        xin = self._alg[self._quantity][Pin]
        xout = self._alg[self._quantity][Pn]
        flow = g*s*(xin - xout)
        
        #Calculate effective property for given algorithm
        if self._alg.__class__.__name__ == 'FickianDiffusion':
            if 'pore.molar_density' in self._fluid.props(mode='scalars'):
                D = sp.sum(flow)*L/A/self._fluid['pore.molar_density']/sp.absolute(BCs[0]-BCs[1])
                return D
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        