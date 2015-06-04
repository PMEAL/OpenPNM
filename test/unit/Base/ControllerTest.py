import OpenPNM
from OpenPNM.Base import Controller
from OpenPNM.Network import GenericNetwork


class ControllerTest:

    def setup_class(self):
        self.controller = Controller()
        self.net = OpenPNM.Network.Cubic(shape=[10,10,10],'test_net')
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=self.net.Ps,
                                                 throats=self.net.Ts)

    def test_str(self):
        actual_string = self.controller.__str__()
        expected_string = \
            '------------------------------------------------------------\n' + \
            'Object:         Name                 (Class)\n' + \
            '------------------------------------------------------------\n' + \
            'Network:        test_net             (Cubic)'

        assert actual_string == expected_string

    def test_get_log_level(self):
        self.controller.loglevel = 50
        assert self.controller.loglevel == 'Log level is currently set to: 50'

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)

    def test_save_and_load(self):
        self.controller.save(self.net.name)
        self.controller.clear()
        assert self.controller == {}
        self.controller.load(self.net.name)
        assert self.net.name in self.controller.keys()
        self.controller.save_simulation(self.net)
        self.controller.purge_object(self.net, mode='complete')
        self.controller.load_simulation(self.net.name)
        self.controller.purge_object(self.net, mode='single')
