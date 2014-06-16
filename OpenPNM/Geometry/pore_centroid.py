r"""
===============================================================================
Submodule -- pore_centroid
===============================================================================

"""
import scipy as sp

def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate the centroid of the pore from the voronoi vertices
    """
    verts = network['pore.vertices']
    value = sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    #network.set_data(prop=propname,pores=geometry.pores(),data=value)
    #network.set_pore_data(prop=propname,data=value)
    network['pore.'+propname]=value