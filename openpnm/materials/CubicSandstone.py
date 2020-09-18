import scipy as sp
import openpnm as op
from openpnm.utils import logging, Project, GenericSettings
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models as mods
logger = logging.getLogger(__name__)


class SandstoneParams():
    whetstone = {'sandstone': 'whetstone',
                 'alpha_t': 0.178e-6,
                 'alpha_b': 0.221e-6,
                 'beta_t': 0.536,
                 'beta_b': 1.18,
                 'bmin_b': 7.1e-6,
                 'bmax_b': 47.5e-6,
                 'bmin_t': 0.55e-6,
                 'bmax_t': 15.0e-6,
                 'lattice_constant': 37.0e-6}
    bandera = {'sandstone': 'bandera',
               'alpha_t': 0.329e-6,
               'alpha_b': 5.097e-6,
               'beta_t': 0.536,
               'beta_b': 1.18,
               'bmin_b': 14.1e-6,
               'bmax_b': 45.8e-6,
               'bmin_t': 1.84e-6,
               'bmax_t': 32.0e-6,
               'lattice_constant': 73.3e-6}
    torpedo = {'sandstone': 'torpedo',
               'alpha_t': 0.755e-6,
               'alpha_b': 2.194e-6,
               'beta_t': 0.536,
               'beta_b': 1.18,
               'bmin_b': 10.6e-6,
               'bmax_b': 45.8e-6,
               'bmin_t': 2.6e-6,
               'bmax_t': 10.0e-6,
               'lattice_constant': 54.6e-6}
    clear_creek = {'sandstone': 'clear_creek',
                   'alpha_t': 0.780e-6,
                   'alpha_b': 3.865e-6,
                   'beta_t': 0.536,
                   'beta_b': 1.18,
                   'bmin_b': 21.1e-6,
                   'bmax_b': 42.3e-6,
                   'bmin_t': 3.321e-6,
                   'bmax_t': 40.0e-6,
                   'lattice_constant': 92.3e-6}
    cottage_grove = {'sandstone': 'cottage_grove',
                     'alpha_t': 0.786e-6,
                     'alpha_b': 3.166e-6,
                     'beta_t': 0.536,
                     'beta_b': 1.18,
                     'bmin_b': 14.1e-6,
                     'bmax_b': 49.3e-6,
                     'bmin_t': 3.0e-6,
                     'bmax_t': 15.0e-6,
                     'lattice_constant': 68.6e-6}
    big_clifty = {'sandstone': 'big_clifty',
                  'alpha_t': 1.137e-6,
                  'alpha_b': 1.870e-6,
                  'beta_t': 0.536,
                  'beta_b': 1.18,
                  'bmin_b': 14.1e-6,
                  'bmax_b': 42.3e-6,
                  'bmin_t': 2.97e-6,
                  'bmax_t': 11.0e-6,
                  'lattice_constant': 71.1e-6}
    st_meinrad = {'sandstone': 'st_meinrad',
                  'alpha_t': 1.269e-6,
                  'alpha_b': 3.903e-6,
                  'beta_t': 0.536,
                  'beta_b': 1.18,
                  'bmin_b': 14.1e-6,
                  'bmax_b': 49.3e-6,
                  'bmin_t': 2.95e-6,
                  'bmax_t': 10.0e-6,
                  'lattice_constant': 69.6e-6}
    noxie_129 = {'sandstone': 'noxie_129',
                 'alpha_t': 0.718e-6,
                 'alpha_b': 2.935e-6,
                 'beta_t': 0.536,
                 'beta_b': 1.18,
                 'bmin_b': 12.3e-6,
                 'bmax_b': 49.3e-6,
                 'bmin_t': 2.7e-6,
                 'bmax_t': 15.0e-6,
                 'lattice_constant': 57.3e-6}
    noxie_47 = {'sandstone': 'noxie_47',
                'alpha_t': 1.842e-6,
                'alpha_b': 4.600e-6,
                'beta_t': 0.536,
                'beta_b': 1.18,
                'bmin_b': 19.4e-6,
                'bmax_b': 52.8e-6,
                'bmin_t': 3.9e-6,
                'bmax_t': 20.0e-6,
                'lattice_constant': 87.0e-6}
    bartlesville = {'sandstone': 'bartlesville',
                    'alpha_t': 2.252e-6,
                    'alpha_b': 4.251e-6,
                    'beta_t': 0.536,
                    'beta_b': 1.18,
                    'bmin_b': 19.0e-6,
                    'bmax_b': 44.1e-6,
                    'bmin_t': 2.69e-6,
                    'bmax_t': 18.5e-6,
                    'lattice_constant': 72.9e-6}
    berea_BE1 = {'sandstone': 'berea_BE1',
                 'alpha_t': 1.697e-6,
                 'alpha_b': 1.803e-6,
                 'beta_t': 0.536,
                 'beta_b': 1.18,
                 'bmin_b': 14.1e-6,
                 'bmax_b': 42.3e-6,
                 'bmin_t': 2.8e-6,
                 'bmax_t': 11.0e-6,
                 'lattice_constant': 73.9e-6}
    berea_108 = {'sandstone': 'berea_108',
                 'alpha_t': 1.904e-6,
                 'alpha_b': 6.081e-6,
                 'beta_t': 0.536,
                 'beta_b': 1.18,
                 'bmin_b': 24.6e-6,
                 'bmax_b': 70.4e-6,
                 'bmin_t': 5.9e-6,
                 'bmax_t': 26.0e-6,
                 'lattice_constant': 125.6e-6}
    boise = {'sandstone': 'boise',
             'alpha_t': 11.050e-6,
             'alpha_b': 8.333e-6,
             'beta_t': 0.536,
             'beta_b': 1.18,
             'bmin_b': 38.7e-6,
             'bmax_b': 78.9e-6,
             'bmin_t': 7.66e-6,
             'bmax_t': 34.0e-6,
             'lattice_constant': 171.2e-6}

class CubicSandstoneSettings(GenericSettings):
    r"""
    The following parameters are used to generate the network and geometry

    The pore and throat sizes are chosen from the Weibull distribution which
    is defined as follows:

    .. math::

        cdf = 1 - exp\bigg[-\bigg(\frac{b_{i} -
              b_{i,min}}{\alpha_{i}}\bigg)^{\beta_{i}}\bigg]

    Values of $\alpha_i$, $beta_i$, and $b_{i,min}$ are unique to each type
    of standstone.  $i=b$ indicates pore body and $i=t$ is throat.

    Parameters
    ----------
    sandstone : str (default = 'berea')
        The type of sandstone to model
    lattice_constant : float
        The lattice spacing used when creating the network
    beta_b and beta_t : float
        Controls the skewness of the Weibull distribution, with smaller
        numbers creating a distribution skewed to the left. In the
        ``scipy.stats.weibull_min`` function this is referred to as ``c``,
        and it is referred to as :math:`\kappa` on the Wikipedia entry
        for the Weibull distribution.  These values are held constant by
        Ioannidis and Chatzis at ``beta_b=1.18`` and ``beta_t=0.536``.
    alpha_b and alpha_t : float
        Controls the scale of the Weibull distribution, with smaller numbers
        creating a distribution shifted to small values.  This is the
        ``scale`` used by the ``scipy.stats.weibull_min`` function
        and is more generally referred to as :math:`\lambda`, for instance
        on the Wikipedia entry for the Weibull distribution.
    bmin_b, bmin_t : float
        Controls the smallest value returned by the Weibull distributions.
    bmax_b, bmax_t : float
        Controls the largest value returned by the Weibull distributions.
    """
    sandstone = 'berea_108'
    lattice_constant = 125.6e-6
    beta_b = 1.18  # The beta values are set to constants for ALL materials
    beta_t = 0.536
    alpha_t = 1.904e-6
    alpha_b = 6.08e-6
    bmin_b = 24.6e-6  # min and max pore sizes are clearly given
    bmax_b = 70.4e-6
    bmin_t = 2.30e-6  # Not sure where this comes from
    bmax_t = 26.4e-6  # Not sure about this, setting to same as pores


class CubicSandstone(Project):
    r"""
    A network representing sandstone using a Cubic lattice.

    This ``material`` is based on the 1993 paper of Ioannidis and Chatzis
    in Chemical Engineering Science entitled "Network modelling of pore
    structure and transport properties of porous media".  The pore and throat
    size parameters are taken from their Table 1 and 5, and are compiled in
    the table below for reference.  To generate the network in the present
    case, these values are looked-up based on the ``sandstone`` parameter
    specified.  They also used a common shape factor for all materials,
    :math:`\beta_b=1.18` and :math:`\beta_t=0.536`.

    Note that all the values listed below are in *um*, while the
    actual values should be given in *m*.

    +---------------+---------+---------+-------+--------+--------+--------+
    | Sandstone     | alpha_t | alpha_b | Lc    | bmin_b | bmax_b | bmin_b |
    +===============+=========+=========+=======+========+========+========+
    | boise         | 13.050  | 8.333   | 171.2 | 38.7   | 73.9   | 1.52   |
    +---------------+---------+---------+-------+--------+--------+--------+
    | bartlesville  |  2.252  | 4.251   |  72.9 | 19.0   | 44.1   | 0.791  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | berea_108     |  1.904  | 6.081   | 125.6 | 24.6   | 70.4   | 0.702  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | berea_BE1     |  1.697  | 8.333   |  73.9 | 14.1   | 42.3   | 0.643  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | noxie_47      |  1.842  | 4.600   |  87.0 | 19.4   | 52.8   | 0.684  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | noxie_129     |  0.718  | 2.935   |  73.3 | 12.3   | 49.3   | 0.296  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | st_meinrad    |  1.269  | 3.903   |  69.6 | 14.1   | 49.3   | 0.503  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | big_clifty    |  1.137  | 1.870   |  71.1 | 14.1   | 42.3   | 0.456  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | cottage_grove |  0.780  | 3.166   |  68.6 | 14.1   | 49.3   | 0.323  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | clear_creek   |  0.780  | 3.865   |  92.3 | 21.1   | 42.3   | 0.321  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | torpedo       |  0.755  | 2.194   |  54.6 | 10.6   | 45.8   | 0.310  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | bandera       |  0.329  | 5.097   |  73.3 | 14.1   | 45.8   | 0.137  |
    +---------------+---------+---------+-------+--------+--------+--------+
    | whetstone     |  0.178  | 0.221   |  37.0 |  7.1   | 47.5   | 0.074  |
    +---------------+---------+---------+-------+--------+--------+--------+

    Parameters
    ----------
    shape : array_like
        The number of pores along each direction of the domain.  All other
        aspects of this model are prescribed by the code.
    sandstone : str
        Options are listed in the table above.  The default is 'berea_108'
    settings : dict
        If ``sandstone`` is given as ``None``, then the user can specify their
        own parameters for the pore and throat size distribution via these
        settings.

    Notes
    -----
    The source code...

    References
    ----------
    [1] ???

    Examples
    --------

    """

    def __init__(self, shape=[10, 10, 10], sandstone='berea_108', settings={},
                 **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(CubicSandstoneSettings())
        self.settings.update(settings)
        if sandstone is not None:
            standstone_settings = getattr(SandstoneParams, sandstone)
            self.settings.update(standstone_settings)
        pn = Cubic(shape=shape, spacing=self.settings['lattice_constant'],
                   connectivity=6, project=self)
# Fetch weibull_min from scipy.stats
        import numpy as np
        import scipy.stats as spst
        wb = spst.weibull_min
# The CDF method tells you the 'probabiity' of finding a value of size 'x'
# or smaller.  By putting the maximum size as 'x', then we're finding the
# maximum seed value that we want to put into the distribution.
        smax_b = wb.cdf(x=self.settings['bmax_b'], c=self.settings['beta_b'],
                        loc=self.settings['bmin_b'],
                        scale=self.settings['alpha_b'])
        smax_t = wb.cdf(x=self.settings['bmax_t'], c=self.settings['beta_t'],
                        loc=self.settings['bmin_t'],
                        scale=self.settings['alpha_t'])
# Start with random numbers between 0 and 1
        geom = GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
        geom['pore.seed'] = np.random.rand(pn.Np)
        geom.add_model(propname='throat.seed',
                     model=op.models.misc.neighbor_lookups.from_neighbor_pores,
                       pore_prop='pore.seed', mode='min')
# Now scale the values by the max found above
        geom['pore.seed'] *= smax_b
        geom['throat.seed'] *= smax_t
#        geom['pore.seed'] = sp.rand(pn.Np) * 0.5 + 0.2
#        geom.add_model(propname='throat.seed',
#                       model=mods.misc.neighbor_lookups.from_neighbor_pores,
#                       pore_prop='pore.seed', mode='min')
        geom.add_model(propname='pore.radius',
                       model=mods.geometry.pore_size.weibull,
                       shape=self.settings['beta_b'],
                       loc=self.settings['bmin_b'],
                       scale=self.settings['alpha_b'],
                       seeds='pore.seed')
        geom.add_model(propname='throat.radius',
                       model=mods.geometry.throat_size.weibull,
                       shape=self.settings['beta_t'],
                       loc=self.settings['bmin_t'],
                       scale=self.settings['alpha_t'],
                       seeds='throat.seed')
        # geom['throat.radius'] = geom['throat.radius']/
        # geom['throat.radius'].max()*45e-6
        # geom['pore.radius'] = geom['pore.radius']/
        # geom['pore.radius'].max()*self.settings['bmax_b']
        # All pores in this model are of square x-section
        # All throats are of slit shape x-section
        geom['pore.size_z'] = geom['pore.radius']*2
        geom['pore.size_x'] = geom['pore.radius']*2
        geom['pore.size_y'] = geom['pore.radius']*3
        geom['throat.size'] = geom['throat.radius']*2

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
        Lc = self.settings['lattice_constant']
        Ts = pn.throats('dir_x')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_x'][Ts] = Lc - geom['pore.size_x'][conns[Ts, 0]]/2\
            - geom['pore.size_x'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with y-directional throats
        Ts = pn.throats('dir_y')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_y'][Ts] = Lc - geom['pore.size_y'][conns[Ts, 0]]/2\
            - geom['pore.size_y'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with z-directional throats
        Ts = pn.throats('dir_z')
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_z'][Ts] = Lc - geom['pore.size_z'][conns[Ts, 0]]/2\
            - geom['pore.size_z'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_z'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_x'][Ts]

        geom.add_model(propname='throat.area',
                       model=mods.misc.basic_math.product,
                       prop1='throat.height', prop2='throat.width')
        geom.add_model(propname='throat.volume',
                       model=mods.misc.basic_math.product,
                       prop1='throat.area', prop2='throat.length')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.basic_math.product,
                       prop1='pore.size_x', prop2='pore.size_y',
                       prop3='pore.size_z')
        import numpy as np
        conns = pn['throat.conns']
        P1, P2 = pn['pore.size_z'][conns.T]
        T = pn['throat.dir_z']
        print(np.all(P1 > T))
        print(np.all(P2 > T))