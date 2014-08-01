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
    This method adds a pore and throat label to the master network (pn) 
    indicating which pores and throats were part of the sub-network (sn).  
    This means that the results of the sub-network can be applied to the 
    correct pores and throats in the main network by calling 
    pn.pores(sn.name), and pn.throats(sn.name).
    
    This also creates a 'pore.map' and 'throat.map' property that lists the 
    correspondence to the subset pores and throats, to the main network.
    
    Examples
    --------
    na
    '''
    newpnm = OpenPNM.Network.GenericNetwork(name=name)
    pores = sp.array(pores,ndmin=1)
    throats = network.find_neighbor_throats(pores=pores,mode='intersection',flatten=True)
    
    #Remap throats on to new pore numbering
    Pmap = sp.zeros_like(network.pores())*sp.nan
    Pmap[pores] = sp.arange(0,sp.shape(pores)[0])
    tpore1 = sp.array(Pmap[network['throat.conns'][throats,0]],dtype=int)
    tpore2 = sp.array(Pmap[network['throat.conns'][throats,1]],dtype=int)
    newpnm['throat.conns'] = sp.vstack((tpore1,tpore2)).T
    newpnm['pore.coords'] = network['pore.coords'][pores]
    newpnm['pore.all'] = network['pore.all'][pores]
    newpnm['throat.all'] = network['throat.all'][throats]
    newpnm['pore.map'] = pores
    newpnm['throat.map'] = throats
    
    #Transfer labels
    labels = network.labels()
    labels.remove('pore.all')
    labels.remove('throat.all')
    for item in labels:
        if item.split('.')[0] == 'pore':
            newpnm[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            newpnm[item] = network[item][throats]
            
    #Convert all Geometry properties into static Network properties
    for geom in network._geometries:
        for prop in geom.props():
            network[prop] = network[prop]
            
    #Transfer props
    props = network.props()
    props.remove('throat.conns')
    props.remove('pore.coords')
    for item in props:
        if item.split('.')[0] == 'pore':
            newpnm[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            newpnm[item] = network[item][throats]
    
    #Append pore and throat mapping to main network as attributes and data
    network['pore.'+newpnm.name] = sp.zeros_like(network['pore.all'],dtype=bool)
    network['pore.'+newpnm.name][pores] = True
    network['throat.'+newpnm.name] = sp.zeros_like(network['throat.all'],dtype=bool)
    network['throat.'+newpnm.name][throats] = True
    
    return newpnm
  
def subset_fluid(fluid,subset_network):
    r'''
    This method takes a fluid from the main network, and converts it to the 
    size and shape of the sub-network.
    
    Parameters
    ----------
    fluid : OpenPNM Fluid Object
        A fluid object that is associated with the main network from which the
        subnetwork was extracted.
        
    subset_network : OpenPNM Network Object
        A subset of a Network object.  The 'pore.map' and 'throat.map' on the
        subset are used to extract data from the correct locations.
        
    Returns
    -------
    newfluid : OpenPNM Fluid Object
        A fluid object with the same shape as the sub-network.  It contains all
        the data of the main fluid, but not the property calculation methods.
        
    Notes
    -----
    This method pulls the data in from the associated Physics objects and adds
    a static copy to the Fluid.
    '''
    newfluid = OpenPNM.Fluids.GenericFluid(network=subset_network)
    for item in fluid.props():
        element = item.split('.')[0]
        newfluid[item] = fluid[item][subset_network[element+'.map']]
    for phys in fluid._physics:
        for item in phys.props():
            element = item.split('.')[0]
            newfluid[item] = phys[item][subset_network[element+'.map']]
    return newfluid


    
    
    