import openpnm.models as mods
import openpnm.models.geometry as gmods
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class SpheresAndCylinders(GenericGeometry):
    r"""


    Parameters
    ----------
    %(GenericGeometry.parameters)s

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

        self.add_model(propname='pore.cross_sectional_area',
                       model=mods.geometry.pore_cross_sectional_area.sphere,
                       pore_diameter='pore.diameter')

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
