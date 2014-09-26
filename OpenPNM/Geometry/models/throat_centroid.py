r"""
===============================================================================
Submodule -- throat_centroid
===============================================================================

"""
import scipy as _sp

def voronoi(geometry,
            **kwargs):
    r"""
    Calculate the centroid of the throat from the voronoi vertices - C.O.M
    """
    verts = geometry['throat.vertices']    
    offset_verts = geometry['throat.offset_vertices']
    value = _sp.ndarray(len(verts),dtype=object)
    for i in range(len(verts)):
        if len(offset_verts[i]) > 2:
            value[i] = _sp.array([offset_verts[i][:,0].mean(),offset_verts[i][:,1].mean(),offset_verts[i][:,2].mean()])
        elif len(verts[i]) > 2:
            value[i] = _sp.array([verts[i][:,0].mean(),verts[i][:,1].mean(),verts[i][:,2].mean()])
        else:
            value[i] = _sp.array([0,0,0])
        
    return value