# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 11:44:11 2014

@author: Stephane
"""
import scipy as sp
from OpenPNM.Visualization import GenericVisualization


class SaveMat(GenericVisualization):
    
    r"""
    SaveMat - Class for writing a mat file to be read by Matlab/Scilab/Octave

    
    Examples
    --------
    Create and store a basic network.

    >>> import OpenPNM
    >>> net = OpenPNM.Network.Cubic(loglevel=20,name='net')
    >>> net.generate(divisions=[20, 20, 20], lattice_spacing=[0.0001],add_boundaries=True)
    >>> vis = OpenPNM.Visualization.SaveMat()
    >>> vis.write(net)

    """
    
    
    def __init__(self,**kwargs):
        r"""
        Initialize
        """

        
        super(SaveMat,self).__init__(**kwargs)
#        self._logger.debug("Execute constructor")
        
    def write(self, network, filename='output.mat', fluids=[]):
        r"""
        Write Network to a VTK file for visualizing in Paraview
    
        Parameters
        ----------
    
        network : OpenPNM Network Object
    
        filename : string ('output.mat')
            Full path to desired file location
            
        fluids : list of fluid objects ([])
            Fluids that have properties we want to write to file
            
    
        """
        pnMatlab = {}        
        new = []
        old = []
        for keys in network.keys():    
            old.append(keys)
            new.append(keys.replace('.','_'))
        
        for i in range(len(network)):        
            pnMatlab[new[i]] = network[old[i]]
                
        
        if len(fluids) != 0:
            for j in range(len(fluids)):
                new = []
                old = []
                
                for keys in fluids[j].keys():    
                    old.append(keys)
                    new.append(fluids[j].name+'_'+keys.replace('.','_'))
                                        
                for i in range(len(fluids[j])):        
                    pnMatlab[new[i]] = fluids[j][old[i]]
            
        sp.io.savemat(file_name=filename,mdict=pnMatlab)
        