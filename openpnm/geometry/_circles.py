from openpnm.models.collections.geometry import circles_and_rectangles
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class CirclesAndRectangles(GenericGeometry):
    r"""

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(circles_and_rectangles)
        self.regenerate_models()
