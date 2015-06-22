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

    def test_save_and_load_workspace(self):
        self.controller.save(self.net.name)
        self.controller.clear()
        assert self.controller == {}
        self.controller.load(self.net.name)
        assert self.net.name in self.controller.keys()

    def test_save_and_load_simulation(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.controller.save_simulation(a)
        assert a in self.controller.values()

    def test_ghost_object(self):
        a = self.controller.ghost_object(self.net)
        # Different objects...
        assert a is not self.net
        # ...but same __dict__ and keys
        assert a.__dict__ == self.net.__dict__
        assert a.keys() == self.net.keys()
        # Ensure an object with same name as a is in Controller dict
        assert a.name in self.controller.keys()
        # But that dictionary key is not a
        assert self.controller[a.name] is not a

    def test_purge_object(self):
        a = OpenPNM.Phases.GenericPhase(network=self.net)
        assert a.name in self.controller.keys()
        assert a in self.controller.values()
        assert a.controller is self.controller
        self.controller.purge_object(a)
        assert a.name not in self.controller.keys()
        assert a not in self.controller.values()
        assert a.controller == {}

    def test_clone_simulation(self):
        a = self.controller.clone_simulation(self.net)
        assert a.name != self.net.name
        assert a in self.controller.values()
        assert a.name in self.controller.keys()

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)
        filelist = os.listdir(".")
        for f in filelist:
            if f.split('.')[-1] in ['pnm', 'net', 'vtp']:
                os.remove(f)


if __name__ == '__main__':
    a = ControllerTest()
    a.setup_class()
    b = a.__class__.__dict__
    for item in b:
        if item.split('_')[0] == 'test':
            print('-'*50)
            print('Testing:'+item)
            b[item](self=a)
    a.teardown_class()
    print('tests complete')
