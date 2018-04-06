import openpnm as op
import scipy as sp
from pathlib import Path
import os


class VTKTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        Ps = [0, 1, 2, 3]
        Ts = self.net.find_neighbor_throats(pores=Ps)
        self.geo_1 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_1['pore.boo'] = 1
        self.geo_1['throat.boo'] = 1
        Ps = [4, 5, 6, 7]
        Ts = self.net.find_neighbor_throats(pores=Ps, mode='intersection')
        self.geo_2 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_2['pore.boo'] = 1
        self.geo_2['throat.boo'] = 1

        self.phase_1 = op.phases.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phases.GenericPhase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2

        self.phys_1 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_1,
                                                geometry=self.geo_1)
        self.phys_1['pore.baz'] = 11
        self.phys_1['throat.baz'] = 11

        self.phys_2 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_1,
                                                geometry=self.geo_2)
        self.phys_2['pore.baz'] = 12
        self.phys_2['throat.baz'] = 12

        self.phys_3 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_2,
                                                geometry=self.geo_1)
        self.phys_3['pore.baz'] = 21
        self.phys_3['throat.baz'] = 21

        self.phys_4 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_2,
                                                geometry=self.geo_2)
        self.phys_4['pore.baz'] = 22
        self.phys_4['throat.baz'] = 22

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()

    def test_save_network(self, tmpdir):
        fname = Path(tmpdir,  'test_save_vtk_1.vtp')
        op.io.VTK.save(network=self.net, filename=fname)
        assert fname.is_file()
        os.remove(fname)

    def test_save_network_and_phases(self, tmpdir):
        fname = Path(tmpdir,  'test_save_vtk_2.vtp')
        op.io.VTK.save(network=self.net, phases=self.phase_1, filename=fname)
        assert fname.is_file()
        os.remove(fname)

    def test_load_no_phases(self, tmpdir):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/VTK-VTP')
        fname = Path(path.resolve(), 'test_save_vtk_1.vtp')
        project = op.io.VTK.load(filename=fname)
        assert len(project) == 1
        net = project.network
        assert net.Np == 8
        assert net.Nt == 12
        assert sp.shape(net['pore.coords']) == (8, 3)
        assert sp.shape(net['throat.conns']) == (12, 2)

    def test_load_with_phases(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/VTK-VTP')
        fname = Path(path.resolve(), 'test_save_vtk_2.vtp')
        project = op.io.VTK.load(filename=fname)
        assert len(project) == 2
        net = project.network
        assert net.Np == 8
        assert net.Nt == 12
        assert sp.shape(net['pore.coords']) == (8, 3)
        assert sp.shape(net['throat.conns']) == (12, 2)
        assert len(project.phases()) == 1


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = VTKTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
