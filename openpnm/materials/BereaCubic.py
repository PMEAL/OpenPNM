from openpnm.utils import logging, Project
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models as mods
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

        geom = GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
        geom['pore.seed'] = sp.rand(pn.Np)
        geom.add_model(propname='throat.seed',
                       model=mods.misc.neighbor_lookups.from_neighbor_pores,
                       pore_prop='pore.seed')

        # Pore throat and pore body characteristic dimensions follow
        # respective Weibull distribution, taking location parameters for
        # Berea 108 sample from table 5.
        geom.add_model(propname = 'pore.size_z',
                       model = mods.geometry.pore_size.weibull,
                       shape = 1.18, loc=6.081e-6, scale = .00004,
                       seeds = 'pore.seed')
        geom.add_model(propname = 'throat.size',
                       model = mods.geometry.throat_size.weibull,
                       shape =0.536, loc = 1.904e-6, scale = .00002,
                       seeds = 'throat.seed')

        # All pores in this model are of square x-section, All throats are of slit shape x-section
        pn['pore.size_x'] = sp.copy(pn['pore.size_z'])
        pn['pore.size_y'] = pn['pore.size_z']*1.5

        # Fetch copies of conns and coords for subsequent size calcs
        conns = pn['throat.conns']
        coords = pn['pore.coords']
        # Create Nt by 2 array of pore coords
        temp = coords[conns]
        temp = sp.absolute(temp[:, 0] - temp[:, 1])

        # Find orientation of each throat and create a label
        pn['throat.dir_x'] = temp[:, 0] > 0
        pn['throat.dir_y'] = temp[:, 1] > 0
        pn['throat.dir_z'] = temp[:, 2] > 0


        # Find width and length of each throat based on it's orientation
        # Start by initializing arrays with 0's
        pn['throat.size_x'] = 0.0
        pn['throat.size_y'] = 0.0
        pn['throat.size_z'] = 0.0
        pn['throat.length'] = 0.0
        pn['throat.width'] = 0.0
        pn['throat.height'] = 0.0

        # Start with x-directional throats
        Ts = pn.throats('dir_x')
        pn['throat.size_z'][Ts] = pn['throat.size'][Ts]
        pn['throat.size_y'][Ts] = pn['throat.size'][Ts]*6
        pn['throat.size_x'][Ts] = Lc - pn['pore.size_x'][conns[Ts, 0]]/2 - pn['pore.size_x'][conns[Ts, 1]]/2
        pn['throat.length'][Ts] = pn['throat.size_x'][Ts]
        pn['throat.width'][Ts] = pn['throat.size_y'][Ts]
        pn['throat.height'][Ts] = pn['throat.size_z'][Ts]

        # Start with y-directional throats
        Ts = pn.throats('dir_y')
        pn['throat.size_z'][Ts] = pn['throat.size'][Ts]
        pn['throat.size_x'][Ts] = pn['throat.size'][Ts]*6
        pn['throat.size_y'][Ts] = Lc - pn['pore.size_y'][conns[Ts, 0]]/2 - pn['pore.size_y'][conns[Ts, 1]]/2
        pn['throat.length'][Ts] = pn['throat.size_y'][Ts]
        pn['throat.width'][Ts] = pn['throat.size_x'][Ts]
        pn['throat.height'][Ts] = pn['throat.size_z'][Ts]

        # Start with z-directional throats
        Ts = pn.throats('dir_z')
        pn['throat.size_x'][Ts] = pn['throat.size'][Ts]
        pn['throat.size_y'][Ts] = pn['throat.size'][Ts]*6
        pn['throat.size_z'][Ts] = Lc - pn['pore.size_z'][conns[Ts, 0]]/2 - pn['pore.size_z'][conns[Ts, 1]]/2
        pn['throat.length'][Ts] = pn['throat.size_z'][Ts]
        pn['throat.width'][Ts] = pn['throat.size_y'][Ts]
        pn['throat.height'][Ts] = pn['throat.size_x'][Ts]

        geom.add_model(propname='throat.cross_sectional_area',
                       model=mods.misc.basic_math.product,
                       prop1='throat.height', prop2='throat.width')
        geom.add_model(propname='throat.volume',
                       model=mods.misc.basic_math.product,
                       prop1='throat.cross_sectional_area', prop2='throat.length')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.basic_math.product,
                       prop1='pore.size_x', prop2='pore.size_y', prop3='pore.size_z')

