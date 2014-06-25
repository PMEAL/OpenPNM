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
    r"""
    """
    def __init__(self,**kwargs):
        r'''
        '''
        super(EffectiveProperty,self).__init__(**kwargs)
        self._logger.info("Construct Algorithm")
        
        
    def setup(self,algorithm,inlets,outlets,conductance=str,quantity=str):
        r'''
        '''
        self._alg = algorithm
        self._conductance = 'throat.'+conductance.split('.')[-1]
        self._quantity = 'pore.'+quantity.split('.')[-1]
        
        #Analyze input and output pores
        #Check for coplanarity
        if self._net.iscoplanar(inlets) == False:
            raise Exception('The specified inlet pores do not define a plane')
        if self._net.iscoplanar(outlets) == False:
            raise Exception('The specified outlet pores do not define a plane')
        #Ensure pores are on a face of domain (only 1 non-self neighbor)
        PnI = self._net.find_neighbor_pores(pores=inlets,mode='not_intersection',excl_self=True)
        if sp.shape(PnI) != sp.shape(inlets):
            raise Exception('The supplied inlet pores have too many neighbors')
            return
        PnO = self._net.find_neighbor_pores(pores=outlets,mode='not_intersection',excl_self=True)
        if sp.shape(PnO) != sp.shape(outlets):
            raise Exception('The supplied outlet pores have too many neighbors')
            return
        
        self._Pin = inlets
        self._Pout = outlets
        
    def run(self):
        r'''
        '''
        Pn = pn.find_neighbor_pores(pores=Pin,excl_self=True)
        Ts = pn.find_connecting_throat(Pin,Pn)
        flow = air['throat.diffusive_conductance'][Ts]*(air['pore.mole_fraction'][Pin] - air['pore.mole_fraction'][Pout])
        
        
        #Fetch area and length of domain
        A = self._net.domain_area(face=self._Pin)
        L = self._net.domain_length(face_1=self._Pin,face_2=self._Pout)
        
        if self._alg.__class__.__name__ == 'FickianDiffusion':
            d1 = self._alg['pore.'+self._quantity][self._Pin][0]
            d2 = self._alg['pore.'+self._quantity][self._Pout][1]
            delta = d1 - d2
#            c = 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        