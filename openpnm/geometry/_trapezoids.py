from openpnm.models.collections.geometry import trapezoids_and_rectangles
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class TrapezoidsAndRectangles(GenericGeometry):
    r"""


    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(trapezoids_and_rectangles)
        self.regenerate_models()
