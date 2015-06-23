import OpenPNM
from os.path import join

class ControllerTest:
    def setup_class(self):
        self.controller = OpenPNM.Base.Controller()
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=self.net.Ps,
                                                 throats=self.net.Ts)

    def test_get_log_level(self):
        self.controller.loglevel = 50
        assert self.controller.loglevel == 'Log level is currently set to: 50'

    def test_save_and_load(self):
        self.controller.save(join(TEMP_DIR, self.net.name))
        self.controller.clear()
        assert self.controller == {}
        self.controller.load(join(TEMP_DIR, self.net.name))
        assert self.net.name in self.controller.keys()

    def test_load_v120_pnm(self):
        temp = self.controller.copy()
        self.controller.clear()
        self.controller.load(join(FIXTURE_DIR, 'test_v120.pnm'))
        a = [
            'Boundary_hy4Ey',
            'FickianDiffusion_LjxxQ',
            'IP_1',
            'OrdinaryPercolation_BI85q',
            'Standard_GIaef',
            'Standard_HmuMH',
            'Toray090_935N3',
            'air',
            'net',
            'water'
        ]
        assert sorted(list(self.controller.keys())) == a
        self.controller.clear()
        self.controller.update(temp)

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)
