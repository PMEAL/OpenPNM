r"""
===============================================================================
Submodule -- throat_length
===============================================================================

"""
import scipy as _sp

def straight(network,
             geometry,
             pore_diameter='pore.diameter',
             **kwargs):
    r"""
    Calculate throat length 
    """
    #Initialize throat_property['length']
    throats = network.throats(geometry.name)
    pore1 = network['throat.conns'][:,0]
    pore2 = network['throat.conns'][:,1]
    C1 = network['pore.coords'][pore1]
    C2 = network['pore.coords'][pore2]
    E = _sp.sqrt(_sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
    D1 = network[pore_diameter][pore1]
    D2 = network[pore_diameter][pore2]
    value = E-(D1+D2)/2
    value = value[throats]
    if _sp.any(value<0):
        geometry._logger.warning('Negative throat lengths are calculated. Arbitrary positive length assigned (1e9 meters)')
        Ts = _sp.where(value<0)[0]
        value[Ts] = 1e-9
    return value
        
def voronoi(network,
            geometry,
            **kwargs):
    r"""
    Calculate the centre to centre distance from centroid of pore1 to centroid of throat to centroid of pore2
    """
    throats = network.throats(geometry.name)
    connections = network['throat.conns'][throats]
    pore1 = connections[:,0]
    pore2 = connections[:,1]
    "This is probably wrong as indices are different"
    pore_centroids = geometry['pore.centroid']
    throat_centroids = geometry['throat.centroid']
    v1 = throat_centroids-pore_centroids[pore1]
    v2 = throat_centroids-pore_centroids[pore2]
    value = _sp.ndarray(len(connections))
    for i in range(len(connections)):
        value[i] = _sp.linalg.norm(v1[i])+_sp.linalg.norm(v2[i])
    return value
