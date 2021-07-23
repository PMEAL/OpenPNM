from openpnm.phases import GenericPhase
from openpnm.models.collections.phase import air
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class Air(GenericPhase):
    r"""

    Parameters
    ----------
    %(GenericPhase.parameters)s


    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(air)
