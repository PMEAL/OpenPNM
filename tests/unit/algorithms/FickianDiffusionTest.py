import openpnm as op
import scipy as sp
import pytest


class GenericTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.phase = op.phases.Air(network=self.net)
        self.phase['pore.mole_fraction'] = 0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1.0

    def test_continuity_BC(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        xyz = self.net["pore.coords"]
        left_interface = (xyz[:, 1] < 5) & (xyz[:, 1] > 4)
        right_interface = (xyz[:, 1] < 6) & (xyz[:, 1] > 5)
        alg.set_continuity_BC(ps1=left_interface, ps2=right_interface, K12=2.0)
        alg.set_value_BC(pores=self.net.pores('left'), values=1.0)
        alg.set_value_BC(pores=self.net.pores('right'), values=0.0)
        alg.run()
        x = alg["pore.mole_fraction"]
        x_LI = sp.unique(sp.around(x[left_interface], decimals=3))
        x_RI = sp.unique(sp.around(x[right_interface], decimals=3))
        assert sp.all(x_LI == 2.0*x_RI)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
