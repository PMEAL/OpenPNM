from openpnm.geometry import GenericGeometry
from openpnm.models.collections.geometry import cones_and_cylinders
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class ConesAndCylinders(GenericGeometry):
    r"""
    blah

    Parameters
    ----------
    %(GenericGeometry.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(cones_and_cylinders)
