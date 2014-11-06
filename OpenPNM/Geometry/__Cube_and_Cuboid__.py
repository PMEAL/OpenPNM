# -*- coding: utf-8 -*-
"""
===============================================================================
Cube_and_Cuboid -- A standard Cubic pore and Cuboic throat model
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry

class Cube_and_Cuboid(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(Cube_and_Cuboid,self).__init__(**kwargs)
        self._generate()

    def _generate(self):
        r'''
        '''
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=self._seed)
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=1.5,
                       psd_loc=14e-6,
                       psd_scale=2e-6)
        self.add_model(propname='pore.area',
                       model=gm.pore_area.cubic)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.cube)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=1.5,
                       tsd_loc=14e-6,
                       tsd_scale=2e-6)
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cuboid)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cuboid)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cuboid)

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

