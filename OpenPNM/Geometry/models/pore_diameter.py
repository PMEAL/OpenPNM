r"""
===============================================================================
pore_diameter -- 
===============================================================================

"""
import scipy as _sp
from numpy.linalg import norm

def sphere(geometry,
           psd_name,psd_shape,psd_loc,psd_scale,
           pore_seed='pore.seed',
           psd_offset=0,
           **kwargs):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    import scipy.stats as spst
    prob_fn = getattr(spst,psd_name)
    P = prob_fn(psd_shape,loc=psd_loc,scale=psd_scale)
    value = P.ppf(geometry[pore_seed])+psd_offset
    return value

def voronoi(geometry,
            pore_volume='pore.volume',
            **kwargs):
    r"""
    Calculate pore diameter from equivalent sphere - volumes must be calculated first
    """
    from scipy.special import cbrt
    pore_vols = geometry[pore_volume]
    value = cbrt(6*pore_vols/_sp.pi)
    return value

def insphere(network,
             geometry,
             **kwargs):
    r"""
    Calculate the diameter representing an inclosed sphere. The maximum is very difficult to caluclate for irregular polygons with more than 4 faces
    so an average distance from the pore centroid to the throat centroid is an approximation
    """
    Np = geometry.num_pores()
    value = _sp.zeros(Np)
    for geom_pore,net_pore in enumerate(geometry["pore.map"]):
        throats = network.find_neighbor_throats(net_pore)
        tcs = network["throat.centroid"][throats]
        pc = network["pore.centroid"][net_pore]
        value[geom_pore]=_sp.mean(norm(tcs-pc))
    return value