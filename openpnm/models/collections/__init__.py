from . import geometry
from . import physics
from . import phase
from . import network
from copy import deepcopy as _deepcopy


class Geometry(object):

    @classmethod
    @property
    def cones_and_cylinders(cls):
        return _deepcopy(geometry.spheres_and_cylinders())
