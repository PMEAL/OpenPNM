"""
module __Cubic__: Generate lattice-like networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""
import OpenPNM

import numpy as np
import scipy as sp
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork

class Cubic(GenericNetwork):
    r"""
    This class contains the methods to create a cubic network with an arbitrary
    domain shape defined by a supplied template.

    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)

    loggername : string
        Overwrite the name of the logger, which defaults to the class name

    Examples
    --------
    This class is a work in progress, examples forthcoming.
    """
    @classmethod
    def from_image(cls, image, *args, **kwargs):
        network = cls(image.shape, *args, **kwargs)
        network['pore.values'] = image.ravel()
        return network
    
    def __init__(self, shape, spacing=None, bbox=None, **kwargs):
        super(Cubic, self).__init__(**kwargs)
        arr = np.atleast_3d(np.empty(shape))

        points = np.array([i for i,v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        if spacing is not None:
            points *= spacing
        elif bbox is not None:
            points *= bbox / self.points.max(axis=0)

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T,H in [
            (I[:,:,:-1], I[:,:,1:]),
            (I[:,:-1], I[:,1:]),
            (I[:-1], I[1:]),
            ]:
            tails.extend(T.flat)
            tails.extend(H.flat)
            heads.extend(H.flat)
            heads.extend(T.flat)
        pairs = np.vstack([tails, heads]).T

        self['pore.coords'] = points
        self['throat.conns'] = pairs

        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)
        self['pore.values'] = arr.ravel()

        x,y,z = self['pore.coords'].T
        self['pore.left'] = x == x.min()
        self['pore.right'] = x == x.max()
        self['pore.bottom'] = y == y.min()
        self['pore.top'] = y == y.max()
        self['pore.back'] = z == z.min()
        self['pore.front'] = z == z.max()

    def add_boundaries(self):
        x,y,z = self['pore.coords'].T
        self['pore.back_face'] = z == z.min()
        self['pore.front_face'] = z == z.max()
        self['pore.bottom_face'] = y == y.min()
        self['pore.top_face'] = y == y.max()
        self['pore.left_face'] = x == x.min()
        self['pore.right_face'] = x == x.max()
        self['pore.boundary'] = x == -1

        t,h = self['throat.conns'].T
        self['throat.back'] = np.ones_like(t, dtype=bool)
        self['throat.back_face'] = np.ones_like(t, dtype=bool)
        self['throat.bottom'] = np.ones_like(t, dtype=bool)
        self['throat.bottom_face'] = np.ones_like(t, dtype=bool)
        self['throat.boundary'] = np.ones_like(t, dtype=bool)
        self['throat.front'] = np.ones_like(t, dtype=bool)
        self['throat.front_face'] = np.ones_like(t, dtype=bool)
        self['throat.left'] = np.ones_like(t, dtype=bool)
        self['throat.left_face'] = np.ones_like(t, dtype=bool)
        self['throat.right'] = np.ones_like(t, dtype=bool)
        self['throat.right_face'] = np.ones_like(t, dtype=bool)
        self['throat.top'] = np.ones_like(t, dtype=bool)
        self['throat.top_face'] = np.ones_like(t, dtype=bool)

    def asarray(self, values):
        points = self['pore.coords']
        spacing = map(np.diff, map(np.unique, points.T))
        min_spacing = [min(a) if len(a) else 1.0 for a in spacing]
        points = (points / min_spacing).astype(int)
        bbox = points.max(axis=0) - points.min(axis=0)
        bbox = (bbox / min_spacing + 1).astype(int)
        actual_indexes = np.ravel_multi_index(points.T, bbox)
        print bbox
        array = np.zeros(bbox)
        array.flat[actual_indexes] = values.ravel()
        return array.squeeze()