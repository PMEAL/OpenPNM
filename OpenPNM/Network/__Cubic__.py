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
    def empty(cls, dims, *a, **kw):
        arr = np.zeros(dims)
        return cls(arr, *a, **kw)

    def __init__(self, ndarray, spacing=None, **kwargs):
        super(Cubic, self).__init__(**kwargs)
        ndarray = np.atleast_3d(ndarray)

        points_rel = np.array(
            [idx for idx,val in np.ndenumerate(ndarray)],
            dtype=float)
        self['pore.coords'] = points_rel

        I = np.arange(ndarray.size).reshape(ndarray.shape)
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
        self['throat.conns'] = np.vstack([tails, heads]).T

        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)
        self['pore.values'] = ndarray.ravel()

        x,y,z = self['pore.coords'].T
        self['pore.back'] = x == x.min()
        self['pore.bottom'] = x == x.min()
        self['pore.front'] = x == x.min()
        self['pore.left'] = x == x.min()
        self['pore.right'] = x == x.min()
        self['pore.top'] = x == x.min()

    def add_boundaries(self):
        x,y,z = self['pore.coords'].T
        self['pore.back_face'] = x == x.min()
        self['pore.bottom_face'] = x == x.min()
        self['pore.boundary'] = x == x.min()
        self['pore.front_face'] = x == x.min()
        self['pore.left_face'] = x == x.min()
        self['pore.right_face'] = x == x.min()
        self['pore.top_face'] = x == x.min()

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

    def asarray(self, values=None):
        # reconstituted facts about the network
        points = self['pore.coords']
        x,y,z = points.T
        span = [(d.max()-d.min() or 1) for d in [x,y,z]]
        res = [len(set(d)) for d in [x,y,z]]

        _ndarray = np.zeros(res)
        rel_coords = np.true_divide(points, span)*(np.subtract(res,1))
        rel_coords = np.rint(rel_coords).astype(int) # absolutely bizarre bug

        actual_indexes = np.ravel_multi_index(rel_coords.T, res)
        if values==None:
            values = self.pores()
        _ndarray.flat[actual_indexes] = values.ravel()
        return _ndarray
