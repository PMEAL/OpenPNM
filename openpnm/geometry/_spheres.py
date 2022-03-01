from openpnm.models.collections.geometry import spheres_and_cylinders
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class SpheresAndCylinders(GenericGeometry):
    r"""
    Pores are treated as spheres and throats as cylinders.

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(spheres_and_cylinders)
        self.regenerate_models()
