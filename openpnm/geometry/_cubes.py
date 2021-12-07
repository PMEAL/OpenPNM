import openpnm.models as mods
import openpnm.models.geometry as gmods
from openpnm.geometry import GenericGeometry


class CubesAndCuboids(GenericGeometry):
    r"""
    Cubes and Cuboids subclass of GenericGeometry. This subclass is
    meant as a basic default geometry to get started quickly.

    Pore diameters are randomly assigned between 0 and the largest sphere
    that does not overlap with it's nearest neighbor.

    Throat diameters are half the diameter of the smaller of it's two
    neighboring pores.

    Parameters
    ----------
    network : GenericNetwork
        The network with which this Geometry should be associated
    pores : array_like
        The pores in the domain where this Geometry applies
    throats : array_like
        The throats in the domain where this Geometry applies
    name : str
        The name of the object, which is also used as the label where this
        geometry is defined.

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
                       props=['pore.max_size', 'pore.seed'])

        self.add_model(propname='pore.volume',
                       model=mods.geometry.pore_volume.cube,
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
                       model=mods.geometry.throat_length.cubes_and_cuboids,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.cross_sectional_area',
                       model=mods.geometry.throat_cross_sectional_area.cuboid,
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.volume',
                       model=mods.geometry.throat_volume.cuboid,
                       throat_diameter='throat.diameter',
                       throat_length='throat.length')

        self.add_model(propname='throat.diffusive_size_factors',
                       model=gmods.diffusive_size_factors.cubes_and_cuboids,
                       pore_diameter="pore.diameter",
                       throat_diameter="throat.diameter")

        self.add_model(propname='throat.hydraulic_size_factors',
                       model=gmods.hydraulic_size_factors.cubes_and_cuboids,
                       pore_diameter="pore.diameter",
                       throat_diameter="throat.diameter")
