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
        if ndims == 1:
            super().__init__(shape=[3, 1, 1], spacing=1e-4)
        elif ndims == 2:
            super().__init__(shape=[3, 3, 1], spacing=1e-4)
        elif ndims == 3:
            super().__init__(shape=[3, 3, 3], spacing=1e-4)
        else:
            raise Exception('ndims must 1, 2 or 3')
        self.add_model_collection(spheres_and_cylinders)
        self.regenerate_models()
