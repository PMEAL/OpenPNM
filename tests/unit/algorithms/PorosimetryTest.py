import openpnm as op
import scipy as sp
import pytest
mgr = op.Workspace()


class PorosimetryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=0.0005)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.water = op.phases.Water(network=self.net)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.water,
                                              geometry=self.geo)
        mod = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.entry_pressure',
                            model=mod)

    def test_(self):
        pass


if __name__ == '__main__':

    t = PorosimetryTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
