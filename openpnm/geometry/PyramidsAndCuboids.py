from openpnm.geometry import GenericGeometry
from openpnm.models.collections.geometry import pyramids_and_cuboids
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class PyramidsAndCuboids(GenericGeometry):
    r"""
    blah

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(pyramids_and_cuboids)
