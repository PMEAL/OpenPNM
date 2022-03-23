import logging
from openpnm.physics import GenericPhysics
from openpnm.models.collections.physics import standard
from openpnm.utils import Docorator


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class Standard(GenericPhysics):
    r"""
    Generic class to generate Physics objects

    Parameters
    ----------
    %(GenericPhysics.parameters)s
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(standard)
        self.regenerate_models()
