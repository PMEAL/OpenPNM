
"""
module MultiPhase
===============================================================================

"""

import OpenPNM
import scipy as sp


def Conduit_Filled_State_Calculator(network):
    r"""
    ---
    """
    if 'Pc_invaded' in network.pore_properties:
        network.pore_properties['swp'] = network.pore_properties['Pc_invaded']>p_val
        network.pore_properties['snwp'] = -network.pore_properties['swp']
    elif 'IP_Pseq' in network.pore_properties:
        network.pore_properties['snwp'] = network.pore_properties['IP_Pseq']>p_val
        network.pore_properties['snwp'] = -network.pore_properties['swp']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    network.throat_properties['swp_conduits'] = -(-network.pore_properties['swp'][pores[:,0]]*-network.pore_properties['swp'][pores[:,1]])
    network.throat_properties['snwp_conduits'] = -(-network.pore_properties['snwp'][pores[:,0]]*-network.pore_properties['snwp'][pores[:,1]])    

#    if -sp.in1d(neighborPs,self.Pinvaded).all():
#        self._net.throat_properties['UninvadedConduits'] = 1
#            val_name = 'Pc_invaded'
#
#        elif self.Alg=='IP':
#            Alg_var = self.Psequence
#            val_name = 'IP_Pseq'
#
#        for P_val in Alg_var:
#            self._logger.info("Applying Pressure/Sequence = "+str(P_val))
#            Concentration = self._do_one_inner_iteration(P_val)
#            #Store result of Concentration in each step
#            if P_val!=0:
#                Concentration = np.multiply(Concentration,self._net.pore_properties[val_name]>P_val)
#      
#            self._net.set_pore_property(name="Concentration_at_"+str(P_val),ndarray=Concentration)
#    return 

def Apply_Phase_State_to_Conduit_Conductivity(network):
    r"""
    ---
    """
    C_wet = network.throat_properties['Cdiff']
    C_wet[network.throat_properties['snwp']] = 1e-30
    return(C_wet)
    
#
#def Late_Pore_Filling_Tracking(network,eta=1,swpstar=0.40):
#    r"""
#    ---
#    """
#    r"""
#    Plot the numerical output of the OP algorithm
#    """
#    PcPoints = np.unique(self._net.pore_properties['Pc_invaded'])
#    Snwp = np.zeros_like(PcPoints)
#    Ps = np.where(self._net.pore_properties['type']==0)
#    for i in range(1,np.size(PcPoints)):
#        Pc = PcPoints[i]
#        Snwp[i] = sum((self._net.pore_properties['Pc_invaded'][Ps]<Pc)*(self._net.pore_properties['volume'][Ps]))/sum(self._net.pore_properties['volume'][Ps])
#        psatnwp[i,0] = sum((1-swpstar*((pn.data['ppcstar'][intPs]/Pc)**eta))*(pn.data['pinvpc'][intPs]<Pc)*(pn.data['pvol'][intPs]))/sum(pn.data['pvol'][intPs])
#        psatnwp[i,1] = sum((pn.data['pinvpc'][intPs]<Pc)*(pn.data['pvol'][intPs]))/sum(pn.data['pvol'][intPs])
#    network.throat_properties['swp'] = swpstar*((network.throat_properties['Pc_entry']<P) /network.self._net.throat_properties['invaded'])**eta
#     print self.indent, '- Late Pore Filling: adjusting saturations'
#        if ('pinvpc' not in pn.data):
#            print self.indent, '! Capillary pressure curve must be simulated first'
#            pn.algs.PCcurve(pn)
#        if ('ppcstar' not in pn.data):
#            print pn.algs.late_pore_filling.__name__ + '>>\t! Pore minimum pressure not given, using Washburn equation'
#            pn.calcPressureFromDiameter('pdia')
#        PcPoints = np.unique(pn.data['pinvpc'])
#        for i in [0,1,2,3]:
#            PcPoints = np.append(PcPoints,PcPoints[-1]*1.25)
#        #Calculate the pore saturation at each pressure step
#        intPs = np.where(pn.data['ptype']==0)
#        psatnwp = np.zeros([np.size(PcPoints),2])
#        pn.data['ppcstar'] = np.min([pn.data['ppcstar'],pn.data['pinvpc']],0)
#        for i in range(1,np.size(PcPoints)):
#            Pc = PcPoints[i]
#            psatnwp[i,0] = sum((1-swpstar*((pn.data['ppcstar'][intPs]/Pc)**eta))*(pn.data['pinvpc'][intPs]<Pc)*(pn.data['pvol'][intPs]))/sum(pn.data['pvol'][intPs])
#            psatnwp[i,1] = sum((pn.data['pinvpc'][intPs]<Pc)*(pn.data['pvol'][intPs]))/sum(pn.data['pvol'][intPs])
#        pn.misc['PcPoints'] = PcPoints
#        pn.misc['satnwp'] = psatnwp
#        plt.plot(PcPoints,psatnwp[:,0],'b.-',PcPoints,psatnwp[:,1],'r.-')
#
#
#    return