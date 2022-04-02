from openpnm.models.collections.phase import water
from openpnm.phase import GenericPhase
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class Water(GenericPhase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    %(GenericPhase.parameters)s

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(water)
        self.regenerate_models()
