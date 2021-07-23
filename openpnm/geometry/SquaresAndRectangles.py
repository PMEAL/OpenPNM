from openpnm.geometry import GenericGeometry
from openpnm.geometry.collections import squares_and_rectangles
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class SquaresAndRectangles(GenericGeometry):
    r"""
    blah

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(squares_and_rectangles)
