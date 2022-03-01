from openpnm.models.collections.geometry import cones_and_cylinders
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class ConesAndCylinders(GenericGeometry):
    r"""
    Pores are treated as cones and throats as cylinders for transport
    conductance. Pores are treated as spheres for all other properties such
    as volume and surface area.

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(cones_and_cylinders)
        self.regenerate_models()
