import scipy as _sp

class Mat():
    
    r"""
    Write Network to a Mat file for exporting to Matlab

    Parameters
    ----------
    
    network : OpenPNM Network Object

    filename : string
        Desired file name, defaults to network name if not given
        
    fluids : list of fluid objects ([])
        Fluids that have properties we want to write to file

    """
    
    def __init__(self,network, filename='', fluids=[],**kwargs):
        r"""
        Initialize
        """
        if filename == '':
            filename = network.name+'.mat'
        self._write(network=network**kwargs)
        
    def _write(self, network, filename='output.mat', fluids=[]):

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
        