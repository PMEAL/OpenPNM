r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as sp

def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate the centroid of the pore from the voronoi vertices - C.O.M
    """
    verts = network['throat.offset_verts']
    value = sp.ndarray(len(verts),dtype=object)
    for i,vert in enumerate(verts):
        value[i] = sp.array([vert[:,0].mean(),vert[:,1].mean(),vert[:,2].mean()])
    network['throat.'+propname]=value