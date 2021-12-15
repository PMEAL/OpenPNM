from openpnm.models.collections.geometry import cones_and_cylinders
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class ConesAndCylinders(GenericGeometry):
    r"""


    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.updates(cones_and_cylinders)
        self.regenerate_models()
