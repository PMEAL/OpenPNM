r"""
===============================================================================
Submodule -- pore_topology
===============================================================================

"""
import scipy as _sp

def get_subscripts(network,
                   shape,
                   **kwargs):
    r'''
    Return the 3D subscripts (i,j,k) into the cubic network
    
    '''
    if network.num_pores('internal') != _sp.prod(shape):
        print('Supplied shape does not match Network size, cannot proceed')
    else:
        template = _sp.atleast_3d(_sp.empty(shape))
        a = _sp.indices(_sp.shape(template))
        i = a[0].flatten()
        j = a[1].flatten()
        k = a[2].flatten()
        ind = _sp.vstack((i,j,k)).T
        vals = _sp.ones((network.Np,3))*_sp.nan
        vals[network.pores('internal')] = ind
        return vals
        
def set_spacing(network,
                spacing,
                subscripts='pore.subscript',
                **kwargs):
    r'''
    Calculates the pore coordinates based on supplied spacing and pore 
    subscript values.  The returned coordinates are offset by half of the 
    lattice spacing so that pore are centered in each lattice cell.  
    
    Notes
    -----
    This model is intended for Cubic networks
    '''
    coords = (network[subscripts] + 0.5)*spacing
    return coords
    








