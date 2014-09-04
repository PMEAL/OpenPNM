r"""
===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy as _sp

def cylinder(geometry,
             tsd_name,tsd_shape,tsd_loc,tsd_scale,
             throat_seed='throat.seed',
             tsd_offset=0,
             **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    import scipy.stats as spst
    prob_fn = getattr(spst,tsd_name)
    P = prob_fn(tsd_shape,loc=tsd_loc,scale=tsd_scale)
    value=P.ppf(geometry[throat_seed])+tsd_offset
    return value
    
def voronoi(geometry,
            throat_area='throat.area',
            **kwargs):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    Equivalent circular diameter from voronoi area
    Could do better here and work out minimum diameter from verts
    """
    areas = geometry[throat_area]
    value = 2*_sp.sqrt(areas/_sp.pi)#64 bit sqrt doesn't work!
    return value
