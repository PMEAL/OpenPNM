from openpnm.geometry import GenericGeometry
from openpnm.geometry.collections import trapezoids_and_rectangles
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class TrapezoidsAndRectangles(GenericGeometry):
    r"""
    blah

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(trapezoids_and_rectangles)

