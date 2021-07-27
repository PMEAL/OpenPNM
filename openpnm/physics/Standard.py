from openpnm.physics import GenericPhysics
from openpnm.models.collections.physics import standard
from openpnm.utils import logging, Docorator, GenericSettings
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class Standard(GenericPhysics):
    r"""
    blah

    Parameters
    ----------
    %(GenericPhysics.parameters)s

    Returns
    -------
    %(GenericPhysics.returns)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(standard)
        self.regenerate_models()
