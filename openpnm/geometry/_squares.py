from openpnm.models.collections.geometry import squares_and_rectangles
from openpnm.geometry import GenericGeometry
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class SquaresAndRectangles(GenericGeometry):
    r"""
    2D representation of cubes and cuboids suitable for 2D simulations


    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(squares_and_rectangles)
        self.regenerate_models()
