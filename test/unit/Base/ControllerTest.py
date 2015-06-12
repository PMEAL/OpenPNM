import OpenPNM
import os
from OpenPNM.Base import Controller


class ControllerTest:

    def setup_class(self):
        self.controller = Controller()
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=self.net.Ps,
                                                 throats=self.net.Ts)

    def test_get_log_level(self):
        self.controller.loglevel = 50
        assert self.controller.loglevel == 'Log level is currently set to: 50'

    def test_save_and_load(self):
        self.controller.save(self.net.name)
        self.controller.clear()
        assert self.controller == {}
        self.controller.load(self.net.name)
        assert self.net.name in self.controller.keys()
        os.remove(self.net.name+'.pnm')

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)
