import OpenPNM
import os
from os.path import join


class WorkspaceTest:
    def setup_class(self):
        self.workspace = OpenPNM.Base.Workspace()
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=self.net.Ps,
                                                 throats=self.net.Ts)

    def test_string(self):
        a = self.workspace.__str__()
        assert type(a) is str

    def test_get_log_level(self):
        self.workspace.loglevel = 50
        assert self.workspace.loglevel == 'Log level is currently set to: 50'

    def test_save_and_load(self):
        self.workspace.save(join(TEMP_DIR, 'test_workspace'))
        self.workspace.clear()
        assert self.workspace == {}
        self.workspace.load(join(TEMP_DIR, 'test_workspace'))
        assert self.net.name in self.workspace.keys()

    def test_load_overwrite_existing(self):
        temp = self.workspace.copy()
        self.workspace.save(join(TEMP_DIR, 'test_workspace'))
        self.workspace.load(join(TEMP_DIR, 'test_workspace'))
        flag = [i for i in temp.keys() if i not in self.workspace.keys()]

    def test_save_no_name(self):
        self.workspace.save()

    def test_load_v120_pnm(self):
        temp = self.workspace.copy()
        self.workspace.clear()
        self.workspace.load(join(FIXTURE_DIR, 'test_v120.pnm'))
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
        assert sorted(list(self.workspace.keys())) == a
        self.workspace.clear()
        self.workspace.update(temp)

    def test_save_and_load_simulation(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.workspace.save_simulation(a, join(TEMP_DIR, 'test_simulation'))
        assert a in self.workspace.values()
        self.workspace.purge_object(a, mode='complete')
        assert a not in self.workspace.values()
        self.workspace.load_simulation(join(TEMP_DIR, 'test_simulation'))
        assert a.name in self.workspace.keys()

    def test_save_simulation_no_name(self):
        a = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.workspace.save_simulation(a)
        self.workspace.clear()
        self.workspace.load_simulation(a.name)

    def test_ghost_object(self):
        a = self.workspace.ghost_object(self.net)
        # Different objects...
        assert a is not self.net
        # ...but same __dict__ and keys
        assert a.__dict__ == self.net.__dict__
        assert a.keys() == self.net.keys()
        # Ensure an object with same name as a is in Workspace dict
        assert a.name in self.workspace.keys()
        # But that dictionary key is not a
        assert self.workspace[a.name] is not a

    def test_purge_object_single(self):
        a = OpenPNM.Phases.GenericPhase(network=self.net)
        assert a.name in self.workspace.keys()
        assert a in self.workspace.values()
        assert a.workspace is self.workspace
        self.workspace.purge_object(a)
        assert a.name not in self.workspace.keys()
        assert a not in self.workspace.values()
        assert a.workspace == {}

    def test_purge_object_complete(self):
        net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net)
        geo.set_locations(pores=net.Ps, throats=net.Ts)
        self.workspace.purge_object(geo, mode='complete')
        assert geo.name not in self.workspace.keys()
        assert net.name not in self.workspace.keys()

    def test_clone_simulation(self):
        a = self.workspace.clone_simulation(self.net)
        assert a.name != self.net.name
        assert a in self.workspace.values()
        assert a.name in self.workspace.keys()

    def test_geometries(self):
        a = self.workspace.geometries()
        assert type(a) is list

    def test_physics(self):
        a = self.workspace.physics()
        assert type(a) is list

    def test_phases(self):
        a = self.workspace.phases()
        assert type(a) is list

    def test_export_VTK_and_MAT(self):
        fname = os.path.join(TEMP_DIR, 'test')
        net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        geo = OpenPNM.Geometry.TestGeometry(network=net,
                                            pores=net.Ps,
                                            throats=net.Ts)
        # Test VTK option
        self.workspace.export(network=net,
                               filename=fname,
                               fileformat='VTK')
        assert os.path.isfile(fname+'.vtp')
        os.remove(fname+'.vtp')
        # Test Matlab matfile option
        self.workspace.export(network=net,
                              filename=fname,
                              fileformat='MAT')
        assert os.path.isfile(fname+'.mat')
        os.remove(fname+'.mat')

    def test_export_one_network_none_specified(self):
        fname = os.path.join(TEMP_DIR, 'test')
        self.workspace.clear()
        OpenPNM.Network.Cubic(shape=[3, 3, 3])
        # Test VTK option
        self.workspace.export(filename=fname,
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
            self.workspace.export(filename=fname,
                                  fileformat='VTK')
        except:
            flag = True
        assert flag

    def test_set_get_comment(self):
        comment = 'Testing the function with a unit test'
        self.workspace.comments = comment
        flag = True
        try:
            self.workspace.comments
        except:
            flag = False
        assert flag

    def teardown_class(self):
        del(self.workspace)
        del(self.net)
        del(self.geo)
