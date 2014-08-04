import OpenPNM
import scipy as sp


def subset_network(network,pores,name=None):
    r'''
    Create a new sub-network from a list of pores.
    
    Parameters
    ----------
    pores : array_like
        A list of pores from which to create the new network
    name : string, optional
        The name to apply to the new network object
        
    Returns
    -------
    OpenPNM Object
        Returns a new network object
        
    Notes
    -----
    This creates a 'pore.map' and 'throat.map' property that lists the 
    correspondence to the subset pores and throats, to the main network.
    
    This method pulls in props from Geometry objects and creates a static copy
    in the new sub-network's dictionary.
    
    Examples
    --------
    na
    '''
    
    subnet = OpenPNM.Network.GenericNetwork(name=name)
    subnet._network = [network]  # Attach parent network
    pores = sp.array(pores,ndmin=1)
    throats = network.find_neighbor_throats(pores=pores,mode='intersection',flatten=True)
    
    #Remap throats on to new pore numbering
    Pmap = sp.zeros_like(network.pores())*sp.nan
    Pmap[pores] = sp.arange(0,sp.shape(pores)[0])
    tpore1 = sp.array(Pmap[network['throat.conns'][throats,0]],dtype=int)
    tpore2 = sp.array(Pmap[network['throat.conns'][throats,1]],dtype=int)
    subnet['throat.conns'] = sp.vstack((tpore1,tpore2)).T
    subnet['pore.coords'] = network['pore.coords'][pores]
    subnet['pore.all'] = network['pore.all'][pores]
    subnet['throat.all'] = network['throat.all'][throats]
    
    #create a map from new pores/throats to parents
    subnet['pore.map'] = pores
    subnet['throat.map'] = throats
    
    #Transfer labels
    labels = network.labels()
    labels.remove('pore.all')
    labels.remove('throat.all')
    for item in labels:
        if item.split('.')[0] == 'pore':
            subnet[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            subnet[item] = network[item][throats]
            
    #Transfer props
    props = network.props()
    props.remove('throat.conns')
    props.remove('pore.coords')
    for geom in network._geometries:
        props.extend(geom.props())
    for item in props:
        if item.split('.')[0] == 'pore':
            subnet[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            subnet[item] = network[item][throats]
    
    return subnet
  
def subset_fluid(fluid,subnet,name=None):
    r'''
    This method takes a fluid from the main network, and converts it to the 
    size and shape of the sub-network.
    
    Parameters
    ----------
    fluid : OpenPNM Fluid Object
        A fluid object that is associated with the main network from which the
        subnetwork was extracted.
        
    subnet : OpenPNM Network Object
        A subset of a Network object.  The 'pore.map' and 'throat.map' on the
        subset are used to extract data from the correct locations.
        
    Returns
    -------
    subfluid : OpenPNM Fluid Object
        A fluid object with the same shape as the sub-network.  It contains all
        the data of the main fluid, but not the property calculation methods.
        
    Notes
    -----
    This method pulls the data in from the associated Physics objects and adds
    a static copy to the new sub-fluid's dictionary.
    '''
    
    subfluid = OpenPNM.Fluids.GenericFluid(network=subnet,name=name)
    pores = subnet['pore.map']
    throats = subnet['throat.map']
    
    #Transfer labels
    labels = fluid.labels()
    labels.remove('pore.all')
    labels.remove('throat.all')
    for item in labels:
        if item.split('.')[0] == 'pore':
            subfluid[item] = fluid[item][pores]
        if item.split('.')[0] == 'throat':
            subfluid[item] = fluid[item][throats]
            
    #Transfer props
    props = fluid.props()
    for phys in fluid._physics:
        props.extend(phys.props())
    for item in props:
        if item.split('.')[0] == 'pore':
            subfluid[item] = fluid[item][pores]
        if item.split('.')[0] == 'throat':
            subfluid[item] = fluid[item][throats]

    return subfluid
    
    
    