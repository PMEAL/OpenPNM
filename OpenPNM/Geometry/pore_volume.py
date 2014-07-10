r"""
===============================================================================
Submodule -- pore_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,pores=geometry.pores(),data=value)

def sphere(geometry,
           network,
           propname,
           diameter='diameter',
           **params):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    value=sp.pi/6*network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def cube(geometry,
         network,
         propname,
         diameter='diameter',
         **params):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    value=network.get_data(prop=diameter,pores=geometry.pores())**3
    network.set_data(prop=propname,pores=geometry.pores(),data=value)
    
def voronoi(geometry,
            network,
            propname,
            **params):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats
    """
    conns = network.get_throat_data(prop='conns')
    verts = network.get_throat_data(prop='offset_verts') 
    num_pores = network.num_pores()
    value = sp.ndarray(num_pores,dtype=object)
    for my_pore in range(num_pores):
        throat_vert_list = []
        num_connections = 0
        for idx,check_pores in enumerate(conns):
            if (check_pores[0] == my_pore) or (check_pores[1] == my_pore):
                num_connections +=1
                for vertex in range(len(verts[idx])):
                    throat_vert_list.append(verts[idx][vertex])
        if num_connections > 1:
            throat_array=sp.asarray(throat_vert_list)
            value[my_pore]= geometry._get_hull_volume(throat_array)
        else:
            value[my_pore]=0.0

    network.set_data(prop=propname,pores=geometry.pores(),data=value)