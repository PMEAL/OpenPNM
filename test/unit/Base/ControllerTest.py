from OpenPNM.Base import Controller
from OpenPNM.Network import GenericNetwork


class ControllerTest:

    def setup_class(self):
        self.controller = Controller()
        self.pn = GenericNetwork(name='test_net')

    def test_str(self):
        actual_string = self.controller.__str__()
        expected_string = \
            '------------------------------------------------------------\n' + \
            'Object:         Name                 (Class)\n' + \
            '------------------------------------------------------------\n' + \
            'Network:        test_net             (GenericNetwork)'

        assert actual_string == expected_string

    def test_get_log_level(self):
        self.controller.loglevel = 50
        assert self.controller.loglevel == 'Log level is currently set to: 50'

    def teardown_class(self):
        del(self.controller)
        del(self.pn)
