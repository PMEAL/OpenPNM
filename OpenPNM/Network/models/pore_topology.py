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
    
    Parameters
    ----------
    shape : list
        The (i,j,k) shape of the network in number of pores in each direction
    
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
        
def apply_spacing(network,
                spacing,
                subscripts='pore.subscript',
                **kwargs):
    r'''
    Calculates the pore coordinates based on supplied spacing and pore 
    subscript values.  The returned coordinates are offset by half of the 
    lattice spacing so that pores are centered in each lattice cell.  
    
    Parameters
    ----------
    spacing : float
        The lattice constant, or spacing between pores.  
    
    Notes
    -----
    - This model is intended for Cubic networks
    - At present this only applied a constant spacing to all directions
    '''
    coords = (network[subscripts] + 0.5)*spacing
    return coords
    
def reduce_coordination(network,
                        z,
                        mode='random',
                        **kwargs):
    r'''
    Reduce the coordination number to the specified z value
    
    Parameters
    ----------
    z : int
        The coordination number or number of throats connected a pore
        
    mode : string, optional
        Controls the logic used to trim connections.  Options are:
        
        - 'random': (default) Throats will be randomly removed to achieve a coordination of z
        - 'max': All pores will be adjusted to have a maximum coordination of z (not implemented yet)
    
    Returns
    -------
    A label array indicating which throats should be trimmed to achieve desired
    coordination.
    
    Notes
    -----
    Pores with only 1 throat will be ignored in all calculations since these 
    are generally boundary pores.
    
    '''
    #Find protected throats
    T_trim = ~network['throat.all']
    P_temp = network.num_neighbors(network.pores())
    T_keep = network.find_neighbor_throats(pores=P_temp==1)
    if mode == 'random':
        T_nums = network.num_neighbors(network.pores())
        z_ave = _sp.average(T_nums[T_nums>1])
        f_trim = (z_ave - z)/z_ave
        T_trim = _sp.rand(network.Nt)<f_trim
        T_trim = T_trim*(~network.tomask(throats=T_keep))
    if mode == 'max':
        pass
    return T_trim






