from openpnm.network import Cubic
from openpnm.utils import Docorator
from openpnm.models.collections.geometry import spheres_and_cylinders


docstr = Docorator()


__all__ = ['Cubic']


@docstr.dedent
class Demo(Cubic):
    r"""

    Parameters
    ----------
    %(GenericNetwork.parameters)s

    """

    def __init__(self, ndims=2, **kwargs):
        if 'shape' in kwargs.keys():
            super().__init__(**kwargs)
        elif ndims == 1:
            super().__init__(shape=[3, 1, 1], spacing=1e-4)
        elif ndims == 2:
            super().__init__(shape=[3, 3, 1], spacing=1e-4)
        elif ndims == 3:
            super().__init__(shape=[3, 3, 3], spacing=1e-4)
        self.add_model_collection(spheres_and_cylinders)
        self.regenerate_models()
