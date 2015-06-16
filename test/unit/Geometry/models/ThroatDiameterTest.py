import OpenPNM
from OpenPNM.Geometry import models as gm
import numpy as np


class ThroatDiameterTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[2, 2, 2])
        self.net.add_boundaries()
        Ps = self.net.pores('boundary', mode='not')
        Ts = self.net.throats('boundary', mode='not')
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=Ps,
                                                 throats=Ts)
        BPs = self.net.pores('boundary')
        BTs = self.net.throats('boundary')
        self.boun = OpenPNM.Geometry.Boundary(network=self.net,
                                              pores=BPs,
                                              throats=BTs)

    def test_min_pore(self):
        self.geo.models.add(propname='throat.diameter',
                            model=gm.throat_diameter.minpore)
        self.boun.models.add(propname='throat.diameter',
                             model=gm.throat_diameter.minpore)
        BPs = self.net["throat.conns"][self.net.throats('boundary')][:, 0]
        assert np.sum(self.net["pore.diameter"][BPs] -
                      self.boun["throat.diameter"]) == 0.0
