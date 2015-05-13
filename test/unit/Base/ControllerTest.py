from OpenPNM.Base import Controller

class ControllerTest:

    def setup_class(self):
      self.controller = Controller()

    def test_str(self):
        pass

    def test_get_log_level(self):
      self.controller.loglevel = 50
      assert self.controller.loglevel == 'Log level is currently set to: 50'

    def teardown_class(self):
      del(self.controller)
