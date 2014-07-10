r"""
===============================================================================
Submodule -- pore_diameter
===============================================================================

"""
import scipy as _sp

def sphere(network,pores,psd_name,psd_shape,psd_loc,psd_scale,**kwargs):
    r"""
    Calculate pore diameter from seed values for a spherical pore body
    """
    import scipy.stats as spst
    prob_fn = getattr(spst,psd_name)
    P = prob_fn(psd_shape,loc=psd_loc,scale=psd_scale)
    value = P.ppf(network['pore.seed'][pores])
    return value

def voronoi(network,pores,**kwargs):
    r"""
    Calculate pore diameter from equivalent sphere - volumes must be calculated first
    """
    from scipy.special import cbrt
    pore_vols = network['pore.volume'][pores]
    value = cbrt(6*pore_vols/_sp.pi)
    return value