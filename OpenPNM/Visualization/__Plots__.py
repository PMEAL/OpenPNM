"""
module __GenericVisualization__: Base class to visualize networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Visualization/__init__.py' file.

"""

import OpenPNM
import scipy as sp

from __GenericVisualization__ import GenericVisualization

class Plots(GenericVisualization):
    r"""
    Methods for plotting common information about pore networks
    
    Parameters
    ----------
    
    Examples
    --------
    >>> print 'nothing yet'
    
    .. note:: 
    n/a
    
    """
    
    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(Plots,self).__init__(**kwargs)
        self._logger.debug("Construct Plots Class")
        
    @staticmethod
    def overview(net):
        r"""
        Plot a montage of key network size distribution histograms
        
        Parameters
        ----------
        net : OpenPNM Network Object
            The network for which the graphs are desired
            
        
        """
        import matplotlib.pyplot as plt
        #Print Pore and Throat size distributions
        plt.clf()
        plt.figure(1)
        plt.subplot(2,2,1)
        plt.hist(net.pore_properties['diameter'][net.pore_properties['type']==0],25,facecolor='green')
        plt.xlabel('Pore Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,2)
        net.get_neighbor_pores(1)
        x = sp.zeros(net.get_num_pores())
        for i in range(0,sp.shape(net.adjacency_matrix['lil']['connections'].rows)[0]):
            x[i] = sp.shape(net.adjacency_matrix['lil']['connections'].rows[i])[0]
        plt.hist(x,25,facecolor='yellow')
        plt.xlabel('Coordination Number')
        plt.ylabel('Frequency')
        plt.subplot(2,2,3)
        plt.hist(net.throat_properties['diameter'][net.throat_properties['type']==0],25,facecolor='blue')
        plt.xlabel('Throat Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,4)
        plt.hist(net.throat_properties['length'][net.throat_properties['type']==0],25,facecolor='red')
        plt.xlabel('Throat Length [m]')
        plt.ylabel('Frequency')
        
    @staticmethod
    def Capillary_Pressure_Curve(net):
        r"""
        Plot drainage capillary pressure curve
        
        Parameters
        ----------
        net : OpenPNM Network Object
            The network for which the graphs are desired
            
        """
        try:
            PcPoints = sp.unique(net.pore_properties['Pc_invaded'])
        except:
            raise Exception('Capillary pressure simulation has not been run')
        import matplotlib.pyplot as plt
        PcPoints = sp.unique(net.pore_properties['Pc_invaded'])
        Snwp = sp.zeros_like(PcPoints)
        Ps = sp.where(net.pore_properties['type']==0)
        for i in range(1,sp.size(PcPoints)):
            Pc = PcPoints[i]
            Snwp[i] = sum((net.pore_properties['Pc_invaded'][Ps]<Pc)*(net.pore_properties['volume'][Ps]))/sum(net.pore_properties['volume'][Ps])
        plt.plot(PcPoints,Snwp,'r.-')
        plt.xlabel('Capillary Pressure')
        plt.ylabel('Invading Fluid Saturation')
        plt.show(block = False)
        