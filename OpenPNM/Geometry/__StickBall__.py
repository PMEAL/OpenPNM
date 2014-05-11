"""
module __StickBall__: Subclass of GenericGeometry for a standard 'stick & ball'
geometrical model
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class Stick_and_Ball(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a 
    basic default geometry to get started quickly.  It's main limitation is 
    that the pore and throat sizes are hard-coded to be weibull distributions
    with specified values.  You *can* override these values by assigning new
    methods to the 'pore_diameter' and 'throat_diameter' attributes using 
    the distribution values you wish.  

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this 
        geometry is defined.
        
    

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Stick_and_Ball,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
   
        self.add_property(prop='pore_seed',model='random')
        self.add_property(prop='throat_seed',model='neighbor_min')
        self.add_property(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc=0,scale=0.5)
        self.add_property(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=0,scale=0.5)
        self.add_property(prop='pore_volume',model='sphere')
        self.add_property(prop='throat_length',model='straight')
        self.add_property(prop='throat_volume',model='cylinder')
        self.add_property(prop='throat_vector',model='pore_to_pore')
        self.add_property(prop='throat_area',model='cylinder')
        self.add_property(prop='throat_surface_area',model='cylinder')
        
if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Geometry.Stick_and_Ball(loglevel=10,name='test_geom',locations=[0],network=pn)
    test.regenerate()
