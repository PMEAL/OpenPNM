from openpnm.utils import logging, Project
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models.geometry as gm
logger = logging.getLogger(__name__)
import scipy as sp


class BereaCubic(Project):
    r"""
    A traditional Berea Sandstone on a Cubic lattice

    Berea Sandstone is one of the standard materials used on geoscience
    studies due to it's importance in oil reservoir engineering as well as
    having well defined pore structure.  This class creates a Cubic Network
    with the appropriate lattice spacing and connectivity, then adds a Geometry
    object with the necessary pore-scale models and prescribed parameters.

    Parameters
    ----------
    shape : array_like
        The number of pores along each direction of the domain.  All other
        aspects of this model are prescribed by the code.

    name : string, optional
        The name to give the Project

    Notes
    -----
    The source code for this Material is relatively straight-forward, so is a
    good example starting point for creating custom materials.

    References
    ----------
    [1] ???

    Examples
    --------

    """

    def __init__(self, shape, name=None, **kwargs):
        super().__init__(name=name)

        Lc = 0.0001256
        pn = Cubic(shape=shape, spacing=Lc, connectivity=6, project=self)
        
        geom = GenericGeometry(network = pn, pores = pn.Ps, throats= pn.Ts)
        geom['pore.seed']= sp.rand(pn.Np)
        # To assign static random seed values between 0 and 1 to each pore and throat of the geometry object
        geom['throat.seed']= sp.rand(pn.Nt)
        
        #Pore throat and pore body characteristic dimensions follow respective Weibull distribution-
        #taking location parameters for Berea 108 sample from table 5.
        
        geom.add_model(propname = 'pore.height', 
                       model = gm.pore_size.weibull, 
                       shape =1.18, loc = 6.081e-6, scale = .00004, 
                       seeds = 'pore.seed')
        
        geom.add_model(propname = 'throat.height', 
                       model = gm.throat_size.weibull, 
                       shape =0.536, loc = 1.904e-6, scale = .00002, 
                       seeds = 'throat.seed')
        
        # All pores in this model are of square x-section, All throats are of slit shape x-section
        
        pn['pore.width'] = pn['pore.height']
        pn['pore.diameter'] = pn['pore.height']
        pn['pore.depth'] = pn['pore.height']* 1.5
        pn['throat.diameter'] = pn['throat.height']
        conns = pn['throat.conns']
        # This is Nt X 2 array which contains indices of the pore pair connected at each end of throat
        coords = pn['pore.coords']
        
        temp = coords[conns]
        # This is Nt X 2 array of (x,y,z) coordinates of pores which lie on each end of the throat
        temp = sp.absolute(temp[:, 0] - temp[:, 1])
        
        pn['throat.xdir'] = temp[:, 0] > 0
        pn['throat.ydir'] = temp[:, 1] > 0
        pn['throat.zdir'] = temp[:, 2] > 0
        
        #pn['throat.height'] = 0.0
        pn['throat.width'] = 0.0
        pn['throat.length'] = 0.0 
        
        # Start with x-directional throats
        Ts = pn.throats('xdir')
        #pn['throat.height'][Ts] = sp.rand()*Lc*0.1
        pn['throat.width'][Ts] = pn['throat.height'][Ts]*6
        pn['throat.length'][Ts] = Lc - pn['pore.width'][conns[Ts, 0]]/2 - pn['pore.width'][conns[Ts, 1]]/2
        
        # Start with y-directional throats
        Ts = pn.throats('ydir')
        #pn['throat.height'][Ts] = sp.rand()*Lc*0.1
        pn['throat.width'][Ts] = pn['throat.height'][Ts]*6
        pn['throat.length'][Ts] = Lc - pn['pore.depth'][conns[Ts, 0]]/2 - pn['pore.depth'][conns[Ts, 1]]/2
        
        # Start with z-directional throats
        Ts = pn.throats('zdir')
        pn['throat.height'][Ts] = sp.rand()*Lc*0.1
        pn['throat.width'][Ts] = pn['throat.height'][Ts]*6
        pn['throat.length'][Ts] = Lc - pn['pore.height'][conns[Ts, 0]]/2 - pn['pore.height'][conns[Ts, 1]]/2
        
        tvol = pn['throat.height']*pn['throat.width']*pn['throat.length']
        pvol = pn['pore.height']*pn['pore.width']*pn['pore.depth']
        poro = (pvol.sum()+tvol.sum())/(1000*(Lc**3))
        pn['pore.volume'] = pvol
        pn['throat.volume'] = tvol
        pn['pore.area'] = 100*(Lc*Lc)
        pn['throat.area'] = pn['throat.height'] * pn['throat.width']
        pn['throat.conduit_lengths.pore1']= 0.0001
        pn['throat.conduit_lengths.throat' ]= 0.00001
        pn['throat.conduit_lengths.pore2' ]=0.0001
