import OpenPNM as op
import scipy as sp
import os
import OpenPNM.Utilities.IO as io


class IOTest:
    def setup_class(self):
        self.net = op.Network.Cubic(shape=[3, 3, 3])
        self.geom = op.Geometry.Cube_and_Cuboid(network=self.net,
                                                pores=self.net.Ps,
                                                throats=self.net.Ts)
        self.phase = op.Phases.Air(network=self.net)
        self.physics = op.Physics.Standard(network=self.net,
                                           phase=self.phase,
                                           pores=self.net.Ps,
                                           throats=self.net.Ts)

    def test_save_load_vtk_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_1')
        io.VTK.save(network=self.net, filename=fname)
        assert os.path.isfile(fname+'.vtp')
        net = io.VTK.load(fname+'.vtp')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)

    def test_save_and_load_vtk_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_2')
        io.VTK.save(network=self.net, filename=fname, phases=self.phase)
        assert os.path.isfile(fname+'.vtp')
        net = io.VTK.load(fname+'.vtp')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert [True for item in net.keys() if 'temperature' in item]
        assert [True for item in net.keys() if 'diffusive_conductance' in item]

    def test_save_and_load_csv_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_csv_1')
        io.CSV.save(network=self.net, filename=fname)
        assert os.path.isfile(fname+'.csv')
        net = io.CSV.load(fname+'.csv')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)

    def test_save_and_load_csv_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_csv_2')
        io.CSV.save(network=self.net, filename=fname, phases=self.phase)
        assert os.path.isfile(fname+'.csv')
        net = io.CSV.load(fname+'.csv')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert [True for item in net.keys() if 'temperature' in item]
        assert [True for item in net.keys() if 'diffusive_conductance' in item]

    def test_save_and_load_mat_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_mat_1')
        io.MAT.save(network=self.net, filename=fname)
        assert os.path.isfile(fname+'.mat')
        net = io.MAT.load(fname+'.mat')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)

    def test_save_and_load_mat_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_mat_2')
        io.MAT.save(network=self.net, filename=fname, phases=self.phase)
        assert os.path.isfile(fname+'.mat')
        net = io.MAT.load(fname+'.mat')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert [True for item in net.keys() if 'temperature' in item]
        assert [True for item in net.keys() if 'diffusive_conductance' in item]
