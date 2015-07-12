import OpenPNM
from os.path import join


class ControllerTest:
    def setup_class(self):
        self.controller = OpenPNM.Base.Controller()
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=self.net.Ps,
                                                 throats=self.net.Ts)

    def test_string(self):
        a = self.controller.__str__()
        assert type(a) is str

    def test_get_log_level(self):
        self.controller.loglevel = 50
        assert self.controller.loglevel == 'Log level is currently set to: 50'

    def test_save_and_load(self):
        self.controller.save(join(TEMP_DIR, 'test_workspace'))
        self.controller.clear()
        assert self.controller == {}
        self.controller.load(join(TEMP_DIR, 'test_workspace'))
        assert self.net.name in self.controller.keys()

    def test_save_and_load_no_filename(self):
        self.controller.save(join(TEMP_DIR, ''))
        self.controller.clear()
        assert self.controller == {}
        for file in os.listdir(TEMP_DIR):
            if file.startswith("20"):  # This will start to fail on Jan 1, 2100
                filename = file
        self.controller.load(join(TEMP_DIR, filename))
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

    def test_save_and_load_simulation(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.controller.save_simulation(a, join(TEMP_DIR, 'test_simulation'))
        assert a in self.controller.values()
        self.controller.purge_object(a, mode='complete')
        assert a not in self.controller.values()
        self.controller.load_simulation(join(TEMP_DIR, 'test_simulation'))
        assert a.name in self.controller.keys()

    def test_save_and_load_simulation_no_filename(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.controller.save_simulation(a, join(TEMP_DIR, ''))
        assert a in self.controller.values()
        self.controller.purge_object(a, mode='complete')
        assert a not in self.controller.values()
        self.controller.load_simulation(join(TEMP_DIR, a.name))
        assert a.name in self.controller.keys()

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

    def test_purge_object_single(self):
        a = OpenPNM.Phases.GenericPhase(network=self.net)
        assert a.name in self.controller.keys()
        assert a in self.controller.values()
        assert a.controller is self.controller
        self.controller.purge_object(a)
        assert a.name not in self.controller.keys()
        assert a not in self.controller.values()
        assert a.controller == {}

    def test_purge_object_complete(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net)
        geo.set_locations(pores=net.Ps, throats=net.Ts)
        self.controller.purge_object(geo, mode='complete')
        assert geo.name not in self.controller.keys()
        assert net.name not in self.controller.keys()

    def test_clone_simulation(self):
        a = self.controller.clone_simulation(self.net)
        assert a.name != self.net.name
        assert a in self.controller.values()
        assert a.name in self.controller.keys()

    def test_geometries(self):
        a = self.controller.geometries()
        assert type(a) is list

    def test_physics(self):
        a = self.controller.physics()
        assert type(a) is list

    def test_phases(self):
        a = self.controller.phases()
        assert type(a) is list

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)
