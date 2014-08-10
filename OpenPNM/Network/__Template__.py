"""
module __Template__: Generate cubic networks from domain templates
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

import numpy as np
import scipy as sp
from scipy.misc import imread
from scipy.ndimage.interpolation import zoom
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork

class Template(GenericNetwork):
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
    def from_file(cls, filename, name=None, dmax=100):
        im = imread(filename)
        im = im[:,:,0]
        rf = np.true_divide(dmax, im.shape).clip(0, 1).min()
        im = zoom(im, rf)

        pn = cls(name=name or os.path.basename(filename), loglevel=0)
        pn.generate(im)

        return pn

    @classmethod
    def from_dir(cls, dirname, ext='.tif', dmax=100):
        filenames = list(sorted(os.path.join(dirname, fn) \
                         for fn in os.listdir(dirname) if ext in fn))
        if not filenames:
            raise Exception("No filenames with '{}' extension found.".format(ext))

        images = []
        for fn in filenames:
            im = imread(fn)
            im = im[:,:,0]
            rf = np.true_divide(dmax, im.shape).clip(0, 1).min()
            im = zoom(im, rf)
            images.append(im)
        stack = np.dstack(images)

        pn = cls(name=os.path.basename(dirname), loglevel=0)
        pn.generate(stack)

        return pn

    def __init__(self, **kwargs):
        super(Template,self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__,": ","Execute constructor")
    
    def generate(self, im):
        r'''
        Parameters
        ----------
        im : 3D array-like
            An image obtained by ie: scipy.misc.imresize
        '''

        if len(im.shape) == 2:
            im = im.reshape(im.shape+(1,))

        # network generation stuff
        coords = np.array([idx for idx,val in np.ndenumerate(im)]).astype(float)

        I = np.arange(im.size).reshape(im.shape)
        heads, tails = [], []
        for A, B in [
            (I[:,:,:-1], I[:,:,1:]),
            (I[:,:-1], I[:,1:]),
            (I[:-1], I[1:]),
            ]:
            hs, ts = np.vstack([A.flat, B.flat])
            heads.append(hs)
            tails.append(ts)
        heads = np.hstack(heads)
        tails = np.hstack(tails)

        # insert into sub-structure
        self['pore.coords'] = coords
        self['pore.all'] = np.ones(len(coords)).astype(bool)
        self['throat.conns'] = np.vstack([heads, tails]).T
        self['throat.all'] = np.ones_like(heads).astype(bool)
        self['pore.values'] = im.ravel()

        # some other labels
        x,y,z = coords.T
        for label, locations in [
            ('left',    x==x.min()),
            ('right',   x==x.max()),
            ('top',     y==y.min()),
            ('bottom',  y==y.max()),
            ('front',   z==z.min()),
            ('back',    z==z.max()),
            ]:
            self.set_pore_info(label=label, locations=locations)

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