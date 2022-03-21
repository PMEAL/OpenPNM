from openpnm.models.collections.geometry import pyramids_and_cuboids
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class PyramidsAndCuboids(GenericGeometry):
    r"""
    Pores are treated as 4-sided pyramids and throats as elongated cuboids
    for transport conductance. Pores are treated as spheres for all other
    properties such as volume and surface area.

    Parameters
    ----------
    %(GenericGeometry.parameters)s


    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(pyramids_and_cuboids)
        self.regenerate_models()
