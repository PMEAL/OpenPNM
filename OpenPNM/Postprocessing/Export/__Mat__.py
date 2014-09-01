import scipy as _sp

class Mat():
    
    r"""
    Write Network to a Mat file for exporting to Matlab

    Parameters
    ----------
    
    network : OpenPNM Network Object

    filename : string
        Desired file name, defaults to network name if not given
        
    phases : list of phase objects ([])
        Phases that have properties we want to write to file

    """
    
    def __init__(self,network, filename='', phases=[],**kwargs):
        r"""
        Initialize
        """
        if filename == '':
            filename = network.name+'.mat'
        self._write(network=network**kwargs)
        
    def _write(self, network, filename='output.mat', phases=[]):

        pnMatlab = {}        
        new = []
        old = []
        for keys in network.keys():    
            old.append(keys)
            new.append(keys.replace('.','_'))
        
        for i in range(len(network)):        
            pnMatlab[new[i]] = network[old[i]]
                
        
        if len(phases) != 0:
            for j in range(len(phases)):
                new = []
                old = []
                
                for keys in phases[j].keys():    
                    old.append(keys)
                    new.append(phases[j].name+'_'+keys.replace('.','_'))
                                        
                for i in range(len(phases[j])):        
                    pnMatlab[new[i]] = phases[j][old[i]]
            
        _sp.io.savemat(file_name=filename,mdict=pnMatlab)
        