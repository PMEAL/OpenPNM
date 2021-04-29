import openpnm.models as mods
import openpnm.models.geometry as gmods
from openpnm.geometry import GenericGeometry


class SpheresAndCylinders(GenericGeometry):
    r"""
    Spheres and Cylinders subclass of GenericGeometry. This subclass is
    meant as a basic default geometry to get started quickly.

    Pore diameters are randomly assigned between 0 and the largest sphere
    that does not overlap with it's nearest neighbor.

    Throat diameters are half the diameter of the smaller of it's two
    neighboring pores.

    Parameters
    ----------
    network : GenericNetwork
        The network with which this Geometry should be associated
    project : Project, optional
        Can be supplied instead of a ``network``
    pores : array_like
        The pores in the domain where this Geometry applies
    throats : array_like
        The throats in the domain where this Geometry applies
    name : str
        The name of the object, which is also used as the label where this
        geometry is defined.

    Examples
    --------
    Geometry objects (along with Physics objects) can be applied to a
    subset of pores and/or throats.  This allows for different geometrical
    property models to be applied in different regions. This is illustrated
    in the following code:

    >>> import numpy as np
    >>> import scipy as sp
    >>> import openpnm as op
    >>> from openpnm.geometry import SpheresAndCylinders
    >>> import matplotlib.pyplot as plt
    >>> pn = op.network.CubicDual(shape=[5, 5, 5])
    >>> Ps = pn.pores('primary')
    >>> Ts = pn.throats('primary')
    >>> geo1 = SpheresAndCylinders(network=pn, pores=Ps, throats=Ts)
    >>> Ps = pn.pores('secondary')
    >>> Ts = pn.throats(['secondary', 'interconnect'])
    >>> geo2 = SpheresAndCylinders(network=pn, pores=Ps, throats=Ts)

    Now override the 'pore.diameter' values on the ``geo2`` object:

    >>> geo2.remove_model('pore.diameter')  # Remove model and data
    >>> geo2['pore.diameter'] = np.random.rand(geo2.Np) * 0.05

    Look at the 'pore.diameter' distributions on each object:

    >>> fig = plt.hist(geo1['pore.diameter'], bins=20, alpha=0.5)
    >>> fig = plt.hist(geo2['pore.diameter'], bins=20, alpha=0.5)

    The resulting figure shows that these two Geometry object each have a
    different pore size distribution, with ``geo2`` being much smaller:

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='pore.seed',
                       model=mods.misc.random,
                       element='pore',
                       num_range=[0.2, 0.7],
                       seed=None)

        self.add_model(propname='pore.max_size',
                       model=mods.geometry.pore_size.largest_sphere,
                       iters=10)

        self.add_model(propname='pore.diameter',
                       model=mods.misc.product,
                       prop1='pore.max_size',
                       prop2='pore.seed')

        self.add_model(propname='pore.volume',
                       model=mods.geometry.pore_volume.sphere,
                       pore_diameter='pore.diameter')

        self.add_model(propname='throat.max_size',
                       model=mods.misc.from_neighbor_pores,
                       mode='min',
                       prop='pore.diameter')

        self.add_model(propname='throat.diameter',
                       model=mods.misc.scaled,
                       factor=0.5,
                       prop='throat.max_size')

        self.add_model(propname='throat.length',
                       model=mods.geometry.throat_length.spheres_and_cylinders,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.cross_sectional_area',
                       model=mods.geometry.throat_cross_sectional_area.cylinder,
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.volume',
                       model=mods.geometry.throat_volume.cylinder,
                       throat_diameter='throat.diameter',
                       throat_length='throat.length')

        self.add_model(propname='throat.diffusive_size_factors',
                       model=gmods.diffusive_size_factors.spheres_and_cylinders,
                       pore_diameter="pore.diameter",
                       throat_diameter="throat.diameter")

        self.add_model(propname='throat.hydraulic_size_factors',
                       model=gmods.hydraulic_size_factors.spheres_and_cylinders,
                       pore_diameter="pore.diameter",
                       throat_diameter="throat.diameter")
