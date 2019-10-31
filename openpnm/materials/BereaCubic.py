import scipy as sp
from openpnm.utils import logging, Project
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models as mods
logger = logging.getLogger(__name__)


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

    def __init__(self, shape=[10, 10, 10], Lc=0.0001256, **kwargs):
        super().__init__(**kwargs)
        pn = Cubic(shape=shape, spacing=Lc, connectivity=6, project=self)
        geom = GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
        geom['pore.seed'] = sp.rand(pn.Np)
        geom.add_model(propname='throat.seed',
                       model=mods.misc.neighbor_lookups.from_neighbor_pores,
                       pore_prop='pore.seed')

        # Pore throat and pore body characteristic dimensions follow
        # respective Weibull distribution, taking location parameters for
        # Berea 108 sample from table 5.
        geom.add_model(propname='pore.size_z',
                       model=mods.geometry.pore_size.weibull,
                       shape=1.85, loc=4.081e-5, scale=0.00001,
                       seeds='pore.seed')
        geom.add_model(propname='throat.size',
                       model=mods.geometry.throat_size.weibull,
                       shape=0.9, loc=1.1e-6, scale=0.000006,
                       seeds='throat.seed')

        # All pores in this model are of square x-section
        # All throats are of slit shape x-section
        geom['pore.size_x'] = sp.copy(geom['pore.size_z'])
        geom['pore.size_y'] = geom['pore.size_z']*1.5

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
        geom['throat.size_x'] = 0.0
        geom['throat.size_y'] = 0.0
        geom['throat.size_z'] = 0.0
        geom['throat.length'] = 0.0
        geom['throat.width'] = 0.0
        geom['throat.height'] = 0.0

        # Start with x-directional throats
        Ts = pn.throats('dir_x')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_x'][Ts] = Lc - geom['pore.size_x'][conns[Ts, 0]] \
            - geom['pore.size_x'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with y-directional throats
        Ts = pn.throats('dir_y')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_y'][Ts] = Lc - geom['pore.size_y'][conns[Ts, 0]]/2 \
            - geom['pore.size_y'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with z-directional throats
        Ts = pn.throats('dir_z')
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_z'][Ts] = Lc - geom['pore.size_z'][conns[Ts, 0]]/2 \
            - geom['pore.size_z'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_z'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_x'][Ts]

        geom.add_model(propname='throat.area',
                       model=mods.misc.basic_math.product,
                       prop1='throat.height', prop2='throat.width')
        geom.add_model(propname='throat.volume',
                       model=mods.misc.basic_math.product,
                       prop1='throat.area',
                       prop2='throat.length')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.basic_math.product,
                       prop1='pore.size_x', prop2='pore.size_y',
                       prop3='pore.size_z')
