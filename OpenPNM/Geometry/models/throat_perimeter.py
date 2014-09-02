r"""
===============================================================================
Submodule -- throat_perimeter
===============================================================================

"""
def voronoi(geometry,
            throat_perimeter='throat.perimeter',
            **kwargs):
    r"""
    As the throat perimeter is stored on the network, lookup the indices pertaining to the geom and retrieve 
    """
    network = geometry._net
    tindex = network.throats(geometry.name)
    value = network[throat_perimeter][tindex]
    return value