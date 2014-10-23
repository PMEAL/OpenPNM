r"""
===============================================================================
throat_offset_vertices -- Offeset throat vertices using a fibre radius parameter
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.vertexops as vo

def voronoi(network,
            geometry,
            offset,
              **kwargs):
    r"""
    Offset the throat vertices effectively erroding the facet by the offset distance supplied
    """    
    Nt = geometry.num_throats()
    area = sp.ndarray(Nt)
    perimeter = sp.ndarray(Nt)
    offset_verts = sp.ndarray(Nt,dtype=object)
    offset_error = sp.ndarray(Nt)
    throat_COM = sp.ndarray([Nt,3])
    for i in range(Nt):
        offset_rand = (sp.random.rand(1)-0.5)*offset            
        throat_verts=geometry["throat.vertices"][i]
        throat_normal=geometry["throat.normal"][i]
        area[i],perimeter[i],offset_verts[i],throat_COM[i],offset_error[i] = vo.get_throat_geom(throat_verts,throat_normal,offset)
    
    for i in range(Nt):
        if offset_error[i] > 0 and len(offset_verts[i]) > 0:
            offset_verts[i]=[]
    "Properties that depend on the offset vertices are the area, perimeter and the centroid or COM"
    "To speed things up we could save them all now rather than processing them individually"
    if kwargs["set_dependent"]==True:
        geometry["throat.area"]=area
        geometry["throat.perimeter"]=perimeter
        geometry["throat.centroid"]=throat_COM
    
    return offset_verts