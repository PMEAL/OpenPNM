import OpenPNM
import os
from os.path import join
import pytest


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

    def test_load_overwrite_existing(self):
        temp = self.controller.copy()
        self.controller.save(join(TEMP_DIR, 'test_workspace'))
        self.controller.load(join(TEMP_DIR, 'test_workspace'))
        flag = [i for i in temp.keys() if i not in self.controller.keys()]

    def test_save_no_name(self):
        self.controller.save()

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
        self.controller.clear()
        assert a not in self.controller.values()
        self.controller.load_simulation(join(TEMP_DIR, 'test_simulation'))
        assert a.name in self.controller.keys()

    def test_save_simulation_no_name(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.controller.save_simulation(a)
        self.controller.clear()
        self.controller.load_simulation(a.name)

    def test_load_simulation_duplicate_names(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10], name='foo')
        b = OpenPNM.Geometry.GenericGeometry(network=a, pores=a.Ps,
                                             throats=a.Ts, name='bar')
        self.controller.save_simulation(a)
        self.controller.clear()
        self.controller.load_simulation(a.name)
        # Will fail since a.name is already present
        with pytest.raises(Exception):
            self.controller.load_simulation(a.name)
        # Update a and b with newly loaded objects
        a = self.controller['foo']
        b = self.controller['bar']
        # Change name of a and it will still fail since name of b is present
        a.name = 'boo'  # Changes name in controller but not file
        with pytest.raises(Exception):
            self.controller.load_simulation(a.name)
        # Change name of b and it will finally pass
        b.name = 'baz'
        self.controller.load_simulation('foo')

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

    def test_import_data(self):
        path = FIXTURE_DIR
        fname = os.path.join(path, 'test_load_csv_no_phases.csv')
        pn = self.controller.import_data(filename=fname)
        assert pn.Np == 27
        fname = os.path.join(path, 'test_load_mat_no_phases.mat')
        pn = self.controller.import_data(filename=fname)
        assert pn.Np == 27
        fname = os.path.join(path, 'test_load_vtk_no_phases.vtp')
        pn = self.controller.import_data(filename=fname)
        assert pn.Np == 27

    def test_export_data(self):
        fname = os.path.join(TEMP_DIR, 'test')
        net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        geo = OpenPNM.Geometry.TestGeometry(network=net,
                                            pores=net.Ps,
                                            throats=net.Ts)
        # Test VTK option
        self.controller.export(network=net,
                               filename=fname,
                               fileformat='VTK')
        assert os.path.isfile(fname+'.vtp')
        os.remove(fname+'.vtp')
        # Test Matlab matfile option
        self.controller.export(network=net,
                               filename=fname,
                               fileformat='MAT')
        assert os.path.isfile(fname+'.mat')
        os.remove(fname+'.mat')
        # Test CSV option
        self.controller.export(network=net,
                               filename=fname,
                               fileformat='csv')
        assert os.path.isfile(fname+'.csv')
        os.remove(fname+'.csv')

    def test_export_one_network_none_specified(self):
        fname = os.path.join(TEMP_DIR, 'test')
        self.controller.clear()
        OpenPNM.Network.Cubic(shape=[3, 3, 3])
        # Test VTK option
        self.controller.export(filename=fname,
                               fileformat='VTK')
        assert os.path.isfile(fname+'.vtp')
        os.remove(fname+'.vtp')

    def test_export_many_networks_none_specified(self):
        OpenPNM.Network.Cubic(shape=[3, 3, 3])
        OpenPNM.Network.Cubic(shape=[3, 3, 3])
        fname = os.path.join(TEMP_DIR, 'test')
        # Test VTK option
        flag = False
        try:
            self.controller.export(filename=fname,
                                   fileformat='VTK')
        except:
            flag = True
        assert flag

    def test_set_get_comment(self):
        comment = 'Testing the function with a unit test'
        self.controller.comments = comment
        flag = True
        try:
            self.controller.comments
        except:
            flag = False
        assert flag

    def teardown_class(self):
        del(self.controller)
        del(self.net)
        del(self.geo)
