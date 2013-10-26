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
import matplotlib.pyplot as plt
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

    def generate(self, stats_pores = {'name' : 'weibull_min',
                                     'shape' : 1.5,
                                       'loc' : 6e-6,
                                     'scale' : 2.5e-5},
                     stats_throats = {'name' : 'weibull_min',
                                     'shape' : 1.5,
                                       'loc' : 6e-6,
                                     'scale' : 2.5e-5},
                          **params):
        r"""
        Generate the network
        """
        self._logger.debug("self.generate()")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        self._add_boundaries()
        self._generate_pore_seeds()
        self._generate_throat_seeds()
        self._generate_pore_diameters(stats_pores)
        self._generate_throat_diameters(stats_throats)
        self._calc_pore_volumes()
        self._calc_throat_lengths()
        self._calc_throat_volumes()
        self._logger.debug("\t end of self.generate()")
        return self._net

    def _generate_setup(self,**params):
        r"""
        Perform applicable preliminary checks and calculations required for generation

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("generation_setup: not implemented")

    def _generate_pores(self):
        r"""
        Generate the pores (numbering, types and coordinates)

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("generate_pores: not implemented")

    def _generate_throats(self):
        r"""
        Generate the throats (numbering and types)

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("generate_throats: not implemented")

    def _add_boundaries(self):
        r"""
        Add boundary pores around network (numbering and types), and create necessary throats

        Notes
        -----
        This method is not implemented in the GenericGeometry and must be sub-classed to produce desired network topology.
        """
        self._logger.error("add_boundaries: not implemented")

    def _generate_pore_seeds(self):
        r"""
        Assign random seed values to pores

        Notes
        -----
        To reproduce an identical network it is necessary to set the seed of the random number generator.  This is accomplished using :code:`scipy.random.seed(<value>)`.
        """
        self._logger.info("generate_pore_seeds: Assign each pore a random seed")
        Np = self._net.get_num_pores()
        self._net.pore_properties['seed'] = sp.random.rand(Np)
        #Overwrite boundary pores with 0 so they don't get assigned values (This won't be necessary with the self-protecting dictionary is used)
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_pore_seeds: End of method")

    def _generate_throat_seeds(self):
        r"""
        Assigns random seeds to throats by adopting the smaller of its neighboring pore seed values

        Notes
        -----
        For this step any boundary pores are assumed to have seed value of 1, which forces the throat to adopt the internal pore's seed.
        """
        self._logger.info("generate_throat_seeds: Assign each throat its smaller neighboring pore seed")
        #Temporarily set boundary pore seeds to 1 to simplify throat seed assignment
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 1
        self._net.throat_properties['seed'] = sp.amin(self._net.pore_properties['seed'][self._net.throat_properties['connections']],1)
        #Set boundary pore seeds back to 0
        self._net.pore_properties['seed'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_throat_seeds: End of method")

    def _generate_pore_diameters(self,stats_pores):
        r"""
        Calculate pore diameter from given statisical distribution using the random seed value for each pore

        Notes
        -----
        The stats_pores dictionary contains the requisite information for the desired distribution.  Each distribution in the Scipy stats library takes slightly different parameters, so this dictionary allows flexibility to send the necessary information.

        """
        self._logger.info("generate_pore_diameters: Generate pore diameter from "+stats_pores['name']+" distribution")
        prob_fn = getattr(spst,stats_pores['name'])
        P = prob_fn(stats_pores['shape'],loc=stats_pores['loc'],scale=stats_pores['scale'])
        self._net.pore_properties['diameter'] = P.ppf(self._net.pore_properties['seed'])
        #Set boundadry pores to size 0
        self._net.pore_properties['diameter'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("generate_pore_diameters: End of method")

    def _generate_throat_diameters(self,stats_throats):
        r"""
        Calculate throat diameter from given statisical distribution using the random seed value for each throat

        Notes
        -----
        The stats_throats dictionary contains the requisite information for the desired distribution.  Each distribution in the Scipy stats library takes slightly different parameters, so this dictionary allows flexibility to send the necessary information.
        """
        self._logger.info("generate_throat_diameters: Generate throat diameter from "+stats_throats['name']+" distribution")
        prob_fn = getattr(spst,stats_throats['name'])
        P = prob_fn(stats_throats['shape'],loc=stats_throats['loc'],scale=stats_throats['scale'])
        self._net.throat_properties['diameter'] = P.ppf(self._net.throat_properties['seed'])
        self._logger.debug("generate_throat_diameters: End of method")

    def _calc_pore_volumes(self):
        r"""
        Calculates pore volume from diameter assuming a spherical pore
        """
        self._logger.info("calc_pore_volumes: Setting pore volumes assuming cubic bodies")
        #Set internal pore volumes to 1
        self._net.pore_properties['volume'] = self._net.pore_properties['diameter']**3
        #Set boundary pore volumes to 0
        self._net.pore_properties['volume'][self._net.pore_properties['type']>0] = 0
        self._logger.debug("calc_pore_volumes: End of method")

    def _calc_throat_volumes(self):
        r"""
        Calculates throat volume from diameter and length assuming a cylindrical pore of constant cross-section

        Notes
        -----
        This calculation does not account for the overlap between the cylinder and spherical pore bodies, so volume is slightly over estimated.
        """
        self._logger.info("calc_throat_volumes: Setting throat volumes assuming square cross-section")
        #Set internal pore volumes to 1
        self._net.throat_properties['volume'] = self._net.throat_properties['length']*self._net.throat_properties['diameter']**2
        self._logger.debug("calc_throat_volumes: End of method")

    def _calc_throat_lengths(self):
        r"""
        Determine throat length from distance between pores minus the radius of each pore

        Notes
        -----
        There is a slight geometric inconsistancy here due to the overlap between the cylinder and the sphereical pore bodies.
        """
        self._logger.info("calc_throat_lengths: Determine throat length from distance between pores")
        #Initialize throat_property['length']
        self._net.throat_properties['length'] = sp.zeros_like(self._net.throat_properties['type'])
        C1 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,0]]
        C2 = self._net.pore_properties['coords'][self._net.throat_properties['connections'][:,1]]
        E = sp.sqrt(sp.sum((C1-C2)**2,axis=1))  #Euclidean distance between pores
        D1 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,0]]
        D2 = self._net.pore_properties['diameter'][self._net.throat_properties['connections'][:,1]]
        self._net.throat_properties['length'] = E - (D1 + D2)/2
        #Perform check for unphysical throat lengths
        if sp.sum(self._net.throat_properties['length']<0):
            self._logger.warning("calc_throat_lengths: Some negative throat lengths exist, some pores overlap!")
        self._logger.debug("calc_throat_lengths: End of method")

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
        print 'not implemented yet'


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

