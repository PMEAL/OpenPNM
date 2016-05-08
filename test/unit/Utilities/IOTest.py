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

    def test_load_statoil(self):
        path = os.path.join(FIXTURE_DIR, 'ICL-SandPack(F42A)')
        net = io.Statoil.load(path=path, prefix='F42A')
        assert net.Np == 1246
        assert net.Nt == 2654
        assert sp.shape(net['pore.coords']) == (1246, 3)
        assert sp.shape(net['throat.conns']) == (2654, 2)
        assert 'pore.radius' in net.keys()

        path = os.path.join(FIXTURE_DIR, 'ICL-Sandstone(Berea)')
        net = io.Statoil.load(path=path, prefix='Berea')
        assert net.Np == 6298
        assert net.Nt == 12098
        assert sp.shape(net['pore.coords']) == (6298, 3)
        assert sp.shape(net['throat.conns']) == (12098, 2)
        assert 'pore.radius' in net.keys()
        assert sp.all(net.find_neighbor_pores(pores=1000) == [221, 1214])

    def test_save_load_vtk_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_1')
        io.VTK.save(network=self.net, filename=fname, legacy=True)
        assert os.path.isfile(fname+'.vtp')
        net = io.VTK.load(fname+'.vtp')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert 'pore.'+self.net.name+'_diameter' in net.keys()

    def test_save_and_load_vtk_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_2')
        io.VTK.save(network=self.net,
                    filename=fname,
                    phases=self.phase,
                    legacy=True)
        assert os.path.isfile(fname+'.vtp')
        net = io.VTK.load(fname+'.vtp')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert [True for item in net.keys() if 'temperature' in item]
        assert [True for item in net.keys() if 'diffusive_conductance' in item]

    def test_save_load_vtk_not_legacy_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_1')
        io.VTK.save(network=self.net,
                    filename=fname,
                    phases=self.phase,
                    legacy=False)
        assert os.path.isfile(fname+'.vtp')
        net = io.VTK.load(fname+'.vtp')
        assert net.Np == 27
        assert net.Nt == 54
        assert sp.shape(net['pore.coords']) == (27, 3)
        assert sp.shape(net['throat.conns']) == (54, 2)
        assert 'pore.diameter' in net.keys()
        assert 'pore.diameter'+'|'+self.net.name not in net.keys()
        assert [item for item in net.keys() if '|'+self.phase.name in item]

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

    def test_load_networkx(self):
        fname = os.path.join(FIXTURE_DIR, 'test_load_yaml.yaml')
        net = io.NetworkX.load(filename=fname)
        assert net.Np == 9
        assert net.Nt == 12
        assert sp.shape(net['pore.coords']) == (9, 3)
        assert sp.shape(net['throat.conns']) == (12, 2)
        a = {'pore.area', 'pore.diameter', 'throat.length', 'throat.perimeter'}
        assert a.issubset(net.props())

    def test_load_imorph(self):
        path = os.path.join(FIXTURE_DIR, 'iMorph-Sandstone')
        path += '/'
        net = io.iMorph.load(node_file=path+'throats_cellsThroatsGraph_Nodes.txt',
                             graph_file=path+'throats_cellsThroatsGraph.txt')
        assert net.Np == 1518
        assert net.Nt == 2424
        assert sp.shape(net['pore.coords']) == (1518, 3)
        assert sp.shape(net['throat.conns']) == (2424, 2)
        a = {'pore.volume', 'pore.types', 'throat.volume', 'throat.types'}
        assert a.issubset(net.props())
        a = {'pore.internal','pore.top_boundary','pore.bottom_boundary',
             'pore.front_boundary','pore.back_boundary','pore.left_boundary',
             'pore.right_boundary'}
        assert a.issubset(net.labels())
