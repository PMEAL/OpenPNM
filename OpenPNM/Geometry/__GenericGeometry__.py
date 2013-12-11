#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import scipy.stats as spst
from functools import partial
import numpy as np

class GenericGeometry(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericGeometry - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)


    """

    def __init__(self, **kwargs):

        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
       
    def create(self,**prms):
        r"""
        Create a geometry object using the supplied parameters
        """
        self.pore_properties = {}
        self.throat_properties = {}
        for key, args in prms.items():
            try:
                function = getattr( getattr(OpenPNM.Geometry, key), args['method'] ) # this gets the method from the file
                preloaded_fn = partial(function, geo=self, **args) 
                setattr(self, key, preloaded_fn)
                print ("Loaded {}.".format(key))
            except AttributeError:
                print( "Did not manage to load {}.".format(key) )
        return self
        
    @staticmethod
    def translate_coordinates(net,displacement=[0,0,0]):
        r"""
        Translate pore network coordinates by specified amount

        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied

        displacement : array_like
            A vector containing the amount to translate in each dimension. [0,0,0] yeilds no translation.

        """
        net.pore_properties['coords'] = net.pore_properties['coords'] + displacement

    @staticmethod
    def scale_coordinates(net,scale=[1,1,1]):
        r"""
        Scale pore network coordinates by specified amount

        Parameters
        ----------
        net : OpenPNM Network Object
            The network to which translation should be applied

        scale : array_like
            A vector containing the amount to scale in each dimension.  [1,1,1] yeilds no scaling.

        """
        net.pore_properties['coords'] = net.pore_properties['coords']*scale

    @staticmethod
    def stitch(net1,net2):
        r"""
        Stitch two networks together

        Parameters
        ----------
        net1 : OpenPNM Network Object
            The network that is stiched to

        net2 : OpenPNM Network Object
            The network that is stitched

        """
        print('not implemented yet')


    def _generate_boundaries(self,**params):
        r"""
        This should generate boundary networks using the parameters passed to generate(), but lattice based on geometry.

        Parameters
        ----------
        params = {
        'psd_info'   : {'name'  :, #Each statistical package takes different params, so send as dict
                'shape' :
                'loc'   :
                'scale' : },
        'tsd_info'   : {'name'  :
                'shape' :
                'loc'   :
                'scale' : },
        'btype'                 :
        'lattice_spacing'       :
        'divisions'             :
        }
        """
        self._logger.error("_generate_boundaries: not implemented")

if __name__ == '__main__':
    test=GenericGeometry(loggername="TestGenerator")

