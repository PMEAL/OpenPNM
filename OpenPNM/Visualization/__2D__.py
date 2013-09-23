# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:01:00 2013

@author: Jeff
"""

import OpenPNM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab



class Vis2D():
    r"""   
    
    Vis2D - Class to plot 2D attributes of the network
    
    Parameters
    ----------
                
    Examples
    --------
    
    
    TODO:
    """
    def overview(self, pn):
        #------------------------------------------------------------------------------
        # These print commands need to be moved to a Vis2D class or something
        #------------------------------------------------------------------------------
        #Print Pore and Throat size distributions
        plt.clf()
        plt.figure(1)
        plt.subplot(2,2,1)
        plt.hist(pn.pore_properties['diameter'][pn.pore_properties['type']==0],25,facecolor='green')
        plt.xlabel('Pore Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,2)
        pn.get_neighbor_pores(1)
        x = np.zeros(pn.get_num_pores())
        for i in range(0,np.shape(pn._adjmatrix_lil.rows)[0]):
            x[i] = np.shape(pn._adjmatrix_lil.rows[i])[0]
        plt.hist(x,25,facecolor='yellow')
        plt.xlabel('Coordination Number')
        plt.ylabel('Frequency')
        plt.subplot(2,2,3)
        plt.hist(pn.throat_properties['diameter'][pn.throat_properties['type']==0],25,facecolor='blue')
        plt.xlabel('Throat Diameter [m]')
        plt.ylabel('Frequency')
        plt.subplot(2,2,4)
        plt.hist(pn.throat_properties['length'][pn.throat_properties['type']==0],25,facecolor='red')
        plt.xlabel('Throat Length [m]')
        plt.ylabel('Frequency')
        #Print Capillary Pressure Curve
#        plt.subplot(3,2,(4,6))
#        OP._plot_results()
#        plt.xlabel('Capillary Pressure [Pa]')
#        plt.ylabel('Non-wetting Phase Saturation')
        
