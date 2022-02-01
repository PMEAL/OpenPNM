from openpnm.physics import GenericPhysics
from openpnm.models.collections.physics import basic
from openpnm.utils import Docorator

docstr = Docorator()
__all__ = ["Basic"]


@docstr.dedent
class Basic(GenericPhysics):
    r"""
    Minimal subclass of GenericPhysics for performing diffusion and/or
    flow simulations.

    Parameters
    ----------
    %(GenericPhysics.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(basic)
        self.regenerate_models()
