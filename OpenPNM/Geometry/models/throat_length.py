r"""
===============================================================================
Submodule -- throat_length
===============================================================================

"""
import scipy as _sp

def straight(network,throats,**kwargs):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    pore1 = network['throat.conns'][throats,0]
    pore2 = network['throat.conns'][throats,1]
    C1 = network['pore.coords'][pore1]
    C2 = network['pore.coords'][pore2]
    E = _sp.sqrt(_sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network['pore.diameter'][pore1]
    D2 = network['pore.diameter'][pore2]
    value = E-(D1+D2)/2
    return value
        
def voronoi(network,throats,**kwargs):
    r"""
    Calculate the centre to centre distance from centroid of pore1 to centroid of throat to centroid of pore2
    """
    connections = network['throat.conns'][throats]
    pore1 = connections[:,0]
    pore2 = connections[:,1]
    pore_centroids = network['pore.centroid']
    throat_centroids = network['throat.centroid']
    v1 = throat_centroids-pore_centroids[pore1]
    v2 = throat_centroids-pore_centroids[pore2]
    value = _sp.ndarray(len(connections), dtype=object)
    for i in range(len(connections)):
        value[i] = _sp.linalg.norm(v1[i])+_sp.linalg.norm(v2[i])
    return value
