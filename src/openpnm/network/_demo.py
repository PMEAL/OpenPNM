from openpnm.network import Cubic
from openpnm.utils import Docorator
from openpnm.models.collections.geometry import spheres_and_cylinders


docstr = Docorator()


__all__ = ['Cubic']


@docstr.dedent
class Demo(Cubic):
    r"""
    A shortcut for generating a cubic network with geometrical properties
    already added.

    Parameters
    ----------
    %(Network.parameters)s

    """

    def __init__(self, shape=[3, 3, 1], **kwargs):
        super().__init__(shape=shape, **kwargs)
        self.add_model_collection(spheres_and_cylinders)
        self.regenerate_models()
