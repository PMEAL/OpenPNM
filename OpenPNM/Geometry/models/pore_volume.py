r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as _sp

def sphere(geometry,
           pore_diameter='pore.diameter',
           **kwargs):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    diams = geometry[pore_diameter]
    value=_sp.pi/6*diams**3
    return value
    
def cube(geometry,
         pore_diameter='pore.diameter',
         **kwargs):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    diams = geometry[pore_diameter]
    value = diams**3
    return value
    
def voronoi(network,
            pores,
            **kwargs):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    """
    conns = network['throat.conns']
    verts = network['throat.offset_verts']
    Np = network.num_pores()
    value = _sp.ndarray(Np,dtype=object)
    for my_pore in range(Np):
        throat_vert_list = []
        num_connections = 0
        for idx,check_pores in enumerate(conns):
            if (check_pores[0] == my_pore) or (check_pores[1] == my_pore):
                num_connections +=1
                for vertex in range(len(verts[idx])):
                    throat_vert_list.append(verts[idx][vertex])
        if num_connections > 1:
            throat_array=_sp.asarray(throat_vert_list)
            value[my_pore]= geometry._get_hull_volume(throat_array)
        else:
            value[my_pore]=0.0
    return value[pores]