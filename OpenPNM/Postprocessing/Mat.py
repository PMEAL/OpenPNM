import scipy as _sp

class save():
    
    r"""
    SaveMat - Class for writing a mat file to be read by Matlab/Scilab/Octave


    """
    
    
    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super().__init__(self,**kwargs)
        self._write(**kwargs)
        
    def _write(self, network, filename='output.mat', fluids=[]):
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
            
        _sp.io.savemat(file_name=filename,mdict=pnMatlab)
        