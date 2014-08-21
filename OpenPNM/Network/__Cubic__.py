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

        x,y,z = self['pore.coords'].T
        self['pore.front'] = x == x.min()
        self['pore.back'] = x == x.max()
        self['pore.left'] = y == y.min()
        self['pore.right'] = y == y.max()
        self['pore.bottom'] = z == z.min()
        self['pore.top'] = z == z.max()

    def add_boundaries(self):
        r'''
        This method uses clone_pore to clone the surface pores (labeled 'left'
        , 'right', etc), then shifts them to the periphery of the domain, and
        gives them the label 'right_face', 'left_face', etc.
        '''
        x,y,z = self['pore.coords'].T
        
        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0,0,0]
        offset['back']  = [x.min(),0,0]
        offset['right'] = [0,y.min(),0]
        offset['top']   = [0,0,z.min()]
        
        scale = {}
        scale['front']  = scale['back']  = [0,1,1]
        scale['left']   = scale['right'] = [1,0,1]
        scale['bottom'] = scale['top']   = [1,1,0]
        
        for label in ['front','back','left','right','bottom','top']:
            ps = self.pores(label)
            self.clone(pores=ps,apply_label=['boundary',label+'_face',label])
            #Translate cloned pores
            ind = self.pores(label+'_face')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords

    def asarray(self, values):
        points = self['pore.coords']
        spacing = map(np.diff, map(np.unique, points.T))
        min_spacing = [min(a) if len(a) else 1.0 for a in spacing]
        points = (points / min_spacing).astype(int)
        bbox = points.max(axis=0) - points.min(axis=0)
        bbox = (bbox / min_spacing + 1).astype(int)
        actual_indexes = np.ravel_multi_index(points.T, bbox)
        print(bbox)
        array = np.zeros(bbox)
        array.flat[actual_indexes] = values.ravel()
        return array.squeeze()