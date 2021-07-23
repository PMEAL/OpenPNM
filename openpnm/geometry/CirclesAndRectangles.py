from openpnm.geometry import GenericGeometry
from openpnm.geometry.collections import circles_and_rectangles
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class CirclesAndRectangles(GenericGeometry):
    r"""
    Blah

    Parameters
    ----------
    %(GenericGeometry.parameters)s


    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.update(circles_and_rectangles)
