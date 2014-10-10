import OpenPNM
import scipy as sp


def subset(network,pores,name=None,**kwargs):
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

    Examples
    --------
    na
    '''
    subnet = OpenPNM.Network.GenericNetwork(name=name,**kwargs)
    subnet._net = network  # Attach parent network
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

    #Transfer labels
    labels = network.labels()
    drop_labels = ['pore.all','throat.all']
    [labels.remove(item) for item in drop_labels]
    for item in labels:
        if item.split('.')[0] == 'pore':
            subnet[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            subnet[item] = network[item][throats]

    #Transfer props
    props = network.props()
    drop_props = ['throat.conns','pore.coords']
    [props.remove(item) for item in drop_props]
    for item in props:
        if item.split('.')[0] == 'pore':
            subnet[item] = network[item][pores]
        if item.split('.')[0] == 'throat':
            subnet[item] = network[item][throats]

    #create a map from new pores/throats to parents
    subnet._logger.debug('Adding pore and throat map to subnet')
    subnet.update({'pore.map' : pores})
    subnet.update({'throat.map' : throats})

    return subnet


