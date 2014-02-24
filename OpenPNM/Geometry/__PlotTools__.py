"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import matplotlib.pyplot as plt


class PlotTools(object):
    r"""
    """

    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(PlotTools,self).__init__(**kwargs)
              
    def show_distributions(self,
                           throat_diameter='diameter',
                           pore_diameter='diameter',
                           throat_length='length',                   
                           fig=None):
        r"""
        Plot a montage of network size distribution histograms
      
        Parameters
        ----------
        fig : Matplotlib figure object, optional
          A pre-existing canvas on which to draw the plots
      
        """
        net = self._net
        if not fig:
            fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax1.hist(net.get_pore_data(prop=pore_diameter)[net.get_pore_indices('all')],25,facecolor='green')
        ax1.set_xlabel('Pore Diameter [m]')
        ax1.set_ylabel('Frequency')

        ax2 = fig.add_subplot(222)
        net.find_neighbor_pores(1)
        x = sp.zeros(net.num_pores())
        for i in list(range(0,sp.shape(net.adjacency_matrix['lil']['connections'].rows)[0])):
            x[i] = sp.shape(net.adjacency_matrix['lil']['connections'].rows[i])[0]
        ax2.hist(x,25,facecolor='yellow')
        ax2.set_xlabel('Coordination Number')
        ax2.set_ylabel('Frequency')

        ax3 = fig.add_subplot(223)
        ax3.hist(net.get_throat_data(prop=throat_diameter)[net.get_throat_indices('all')],25,facecolor='blue')
        ax3.set_xlabel('Throat Diameter [m]')
        ax3.set_ylabel('Frequency')
        
        ax4 = fig.add_subplot(224)
        ax4.hist(net.get_throat_data(prop=throat_length)[net.get_throat_indices('all')],25,facecolor='red')
        ax4.set_xlabel('Throat Length [m]')
        ax4.set_ylabel('Frequency')

if __name__ == '__main__':
    pass
