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
from scipy.misc import imread, imresize
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
    >>> img = sp.ones((30,30,30),dtype=int)
    >>> pn = OpenPNM.Network.Template(name='template_1').generate(template=img,lattice_spacing=0.001)
    >>> pn.num_pores()
    27000
    >>> pn.num_throats()
    78300
    """

    def __init__(self, **kwargs):
        super(Template,self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__,": ","Execute constructor")
    
    def generate(self, im, dmax=200, threshold=None):
        r'''
        Parameters
        ----------
        im : 2D or 3D array-like
            A monochrome image obtained by ie: scipy.misc.imresize
        dmax : float
            The maximum possible number of voxels in any dimension.
            Useful for downsizing large networks for preview
        threshold : float
            Images are currently assumed to be dark in areas of interest,
            and clear in void areas. The threshold value determines which
            pores and throats will be pruned out.

            e.g.: with a threshold of 0.5,
                    0 0 0       o-o-o
                    1 1 1   ->     
                    0 0 0       o-o-o
        '''
        assert(len(im.shape) in [2,3])
        im = self._process_image(im, dmax)

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

        if threshold:
            # prune bad
            accessible = I[im < threshold]
            good_heads = np.in1d(heads, accessible)
            good_tails = np.in1d(tails, accessible)
            heads = heads[good_heads & good_tails]
            tails = tails[good_heads & good_tails]

            # every id in tails maps somewhere in accessible
            coords = coords[accessible]
            translate = dict(zip(accessible, np.arange(accessible.size)))
            heads = np.array(map(translate.get, heads))
            tails = np.array(map(translate.get, tails))

        # insert into sub-structure
        self.set_pore_data(prop='coords', data=coords)
        self.set_pore_info(label='all', locations=np.ones(len(coords)).astype(bool))
        self.set_throat_data(prop='connections', data=np.vstack([heads, tails]).T)
        self.set_throat_info(label='all', locations=np.ones_like(heads).astype(bool))

    @staticmethod
    def _process_image(im, dmax):
        rf = np.true_divide(dmax, im.shape).clip(0,1).min()
        if rf != 1:
            im = imresize(im, [int(d*rf) for d in im.shape]) # downsize large
        im = im.astype(float)
        im-= im.min() # normalize
        im/= im.max()
        if len(im.shape)==2:
            im = im.T
            im = im.reshape(im.shape+(1,))
        return im 

    def add_pore_property_from_template(self, template, prop):
        r"""
        Add pore properties based on value stored at each location in the template array

        Parameters
        ----------
        template : array_like
            The template array containing the pore property values at the desired locations

        prop : string
            The name of the pore property being added

        Notes
        -----
        This method can lead to troubles if not executed in the right order.
        For instance, if throat sizes are assigned during the generation stage
        based on neighboring pores sizes, then rewriting pore sizes with this
        method could invalidate the throat sizes.  Ideally, this method should
        be called by the generate_pore_sizes() step during the generate process
        to avoid the issue.  Alternatively, an update_throat_sizes() method
        could be written and called subsequent to calling this method.

        """
        self._logger.info("add_pore_prop_from_template: Add pore properties")
        pore_prop = sp.ravel(template)[self.get_pore_data(prop='voxel_index')]
        self.set_pore_data(prop=prop, data=pore_prop)
        self._logger.debug("add_pore_prop_from_template: End of method")