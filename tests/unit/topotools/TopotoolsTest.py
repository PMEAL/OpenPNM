import numpy as np
import pytest
from numpy.testing import assert_allclose

import openpnm as op
from openpnm import topotools


class TopotoolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_filter_pores_by_z(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        Ps = op.topotools.filter_pores_by_z(network=pn, pores=0, z=1)
        assert_allclose(Ps, [])
        Ps = op.topotools.filter_pores_by_z(network=pn, pores=0, z=2)
        assert_allclose(Ps, [0])
        Ps = op.topotools.filter_pores_by_z(network=pn, pores=[0, 1], z=1)
        assert_allclose(Ps, [])
        Ps = op.topotools.filter_pores_by_z(network=pn, pores=[0, 1], z=2)
        assert_allclose(Ps, [0])
        Ps = op.topotools.filter_pores_by_z(network=pn, pores=pn.Ps, z=4)
        assert_allclose(Ps, [5, 6, 9, 10])

    def test_label_faces(self):
        net = op.network.Cubic(shape=[3, 3, 3], connectivity=6)
        net.clear(mode='labels')
        assert net.labels() == []
        topotools.label_faces(network=net)
        assert net.num_pores('surface') == 26
        assert net.num_pores('left') == 9
        assert net.num_pores('right') == 9
        assert net.num_pores('front') == 9
        assert net.num_pores('back') == 9
        assert net.num_pores('top') == 9
        assert net.num_pores('bottom') == 9

    def test_label_faces_tol(self):
        net = op.network.Cubic(shape=[3, 3, 3], spacing=1, connectivity=6)
        net.clear(mode='labels')
        net['pore.coords'] += np.array([5, 5, 5])
        topotools.label_faces(network=net, tol=0.2)
        assert net.num_pores('surface') == 26
        assert net.num_pores('left') == 9
        assert net.num_pores('right') == 9
        assert net.num_pores('front') == 9
        assert net.num_pores('back') == 9
        assert net.num_pores('top') == 9
        assert net.num_pores('bottom') == 9

    def test_find_surface_pores_default_markers(self):
        from skimage.morphology import ball
        net = op.network.CubicTemplate(template=ball(3), spacing=1)
        net.clear(mode='labels')
        assert net.labels() == []
        topotools.find_surface_pores(network=net)
        assert net.num_pores('surface') == 66

    def test_find_surface_pores_custom_markers_2d(self):
        net = op.network.Cubic(shape=[4, 4, 1], spacing=1)
        net.clear(mode='labels')
        assert net.labels() == []
        markers = [[-1, 2], [2, -1], [2, 5], [5, 2]]
        topotools.find_surface_pores(network=net, markers=markers)
        assert net.num_pores('surface') == 12
        markers = [[-1], [2], [2], [5]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)
        markers = [[-1, 2, 0], [2, -1, 0], [2, 5, 0], [5, 2, 0]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)

    def test_find_surface_pores_custom_markers_3d(self):
        net = op.network.Cubic(shape=[4, 4, 4], spacing=1)
        net.clear(mode='labels')
        assert net.labels() == []
        markers = [[-1, 2, 2], [2, -1, 2], [2, 5, 2], [5, 2, 2]]
        topotools.find_surface_pores(network=net, markers=markers)
        assert net.num_pores('surface') == 48
        markers = [[-1], [2], [2], [5]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)
        markers = [[-1, 2], [2, -1], [2, 5], [5, 2]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)

    def test_template_sphere_shell(self):
        im = topotools.template_sphere_shell(r_outer=4, r_inner=2)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 218
        assert net.Nt == 480
        im = topotools.template_sphere_shell(r_outer=4)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 251
        assert net.Nt == 618

    def test_template_cylinder_annulus(self):
        im = topotools.template_cylinder_annulus(z=10, r_outer=4, r_inner=2)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 320
        assert net.Nt == 688
        im = topotools.template_cylinder_annulus(z=10, r_outer=4)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 450
        assert net.Nt == 1165
        im = topotools.template_cylinder_annulus(z=1, r_outer=4)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 45
        assert net.Nt == 76

    def test_add_boundary_pores(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        topotools.add_boundary_pores(network=net, pores=net.pores('left'),
                                     offset=[0, 1, 0])
        assert net.Np == 150

    def test_clone_pores_mode_parents(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        topotools.clone_pores(network=net, pores=net.pores('left'))
        assert net.Np == 150
        assert net.Nt == 325

    def test_clone_pores_mode_sibings(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        topotools.clone_pores(network=net, pores=net.pores('left'),
                              mode='siblings')
        assert net.Np == 150
        assert net.Nt == 340

    def test_clone_pores_mode_isolated(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        topotools.clone_pores(network=net, pores=net.pores('left'),
                              mode='isolated')
        assert net.Np == 150
        assert net.Nt == 300

    def test_clone_pores_with_labels(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        topotools.clone_pores(network=net, pores=net.pores('left'),
                              labels=['test1', 'test2'])
        assert net.num_pores('test1') == 25
        assert net.num_pores('test2') == 25

    def test_merge_networks(self):
        net1 = op.network.Cubic(shape=[3, 3, 3])
        net2 = op.network.Cubic(shape=[3, 3, 3])
        net1['pore.test1'] = True
        net1['pore.test2'] = 10
        net1['pore.test3'] = np.ones((net1.Np, 3))
        net2['pore.test4'] = True
        net2['pore.test5'] = 10.0
        net2['pore.test6'] = np.ones((net2.Np, 2))
        topotools.merge_networks(network=net1, donor=net2)
        assert np.sum(net1['pore.test1']) == 27
        assert np.all(net1['pore.test3'].shape == (54, 3))
        assert np.sum(net1['pore.test2'][:27]) == 270.0
        assert np.sum(net1['pore.test4'][27:]) == 27
        assert 'pore.test1' not in net2
        assert 'pore.test2' not in net2

    def test_merge_networks_with_active_geometries(self):
        pn = op.network.Cubic(shape=[3, 3, 3], name='net_01')
        pn2 = op.network.Cubic(shape=[3, 3, 3], name='net_02')
        pn2['pore.coords'] += [0, 0, 3]
        pn2['pore.net_02'] = True
        pn2['throat.net_02'] = True
        pn2['pore.test_vals'] = 1.5
        op.topotools.merge_networks(network=pn, donor=pn2)
        assert isinstance(pn['pore.test_vals'], np.ndarray)
        assert ('pore.' + pn2.name) in pn.keys()

    def test_merge_pores(self):
        testnet = op.network.Cubic(shape=[10, 10, 10])
        to_merge = [[0, 1], [998, 999]]
        topotools.merge_pores(testnet, to_merge)
        assert testnet.Np == 998

    def test_connect_pores(self):
        testnet = op.network.Cubic(shape=[10, 10, 10])
        Nt_old = testnet.Nt
        ps1 = [[0, 1], [23, 65]]
        ps2 = [[55], [982, 555]]
        topotools.connect_pores(testnet, pores1=ps1, pores2=ps2)
        am = testnet.create_adjacency_matrix(weights=np.ones(testnet.Nt,
                                                             dtype=int),
                                             fmt='csr')
        conns = testnet['throat.conns']
        assert len(conns) == Nt_old + 6
        assert am[0, 55] == 1
        assert am[1, 55] == 1
        assert am[23, 982] == 1
        assert am[23, 555] == 1
        assert am[65, 982] == 1
        assert am[65, 555] == 1

    def test_trim_pores(self):
        np.random.seed(1)
        pn = op.network.Cubic(shape=[2, 2, 2], spacing=1)
        Ps = pn.pores()[:4]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        pn.set_label(pores=Ps, throats=Ts, label='domain1')
        Ps = pn.pores()[4:]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='union')
        pn.set_label(pores=Ps, throats=Ts, label='domain2')
        pn['pore.random@domain1'] = np.random.random(pn.num_pores('domain1'))
        pn['pore.random@domain2'] = np.random.random(pn.num_pores('domain2'))
        trimmers = pn['pore.random'] < 0.25
        topotools.trim(pn, pores=pn.pores()[trimmers])
        assert ~np.any(pn['pore.random'] < 0.25)

    def test_trim_throats(self):
        np.random.seed(1)
        pn = op.network.Cubic(shape=[2, 2, 2], spacing=1)
        Ps = pn.pores()[:4]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        pn.set_label(pores=Ps, throats=Ts, label='domain1')
        Ps = pn.pores()[4:]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='union')
        pn.set_label(pores=Ps, throats=Ts, label='domain2')
        pn['throat.random@domain1'] = np.random.random(pn.num_throats('domain1'))
        pn['throat.random@domain2'] = np.random.random(pn.num_throats('domain2'))
        trimmers = pn['throat.random'] < 0.25
        topotools.trim(pn, throats=pn.Ts[trimmers])
        assert ~np.any(pn['throat.random'] < 0.25)

    def test_iscoplanar(self):
        # Generate planar points with several parallel vectors at start
        coords = [[0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 2]]
        assert topotools.iscoplanar(coords)
        # NON-planar points, also with parallel vectors
        coords = [[0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 2], [1, 1, 2]]
        assert ~topotools.iscoplanar(coords)
        # Planar points, none parallel
        coords = [[0, 0, 0], [0, 1, 2], [0, 2, 1], [0, 3, 2], [0, 2, 3]]
        assert topotools.iscoplanar(coords)
        # Non-planar points, none parallel
        coords = [[0, 0, 0], [0, 1, 2], [0, 2, 1], [0, 3, 3], [1, 1, 2]]
        assert ~topotools.iscoplanar(coords)

    def test_extend(self):
        pn = op.network.Cubic(shape=[2, 2, 1])
        pn['pore.test_float'] = 1.0
        pn['pore.test_int'] = 1
        pn['pore.test_bool'] = True
        op.topotools.extend(network=pn, coords=[[3, 3, 3], [3, 3, 4]])
        assert np.any(np.isnan(pn['pore.test_float']))
        assert np.any(np.isnan(pn['pore.test_int']))
        assert pn['pore.test_bool'].sum() < pn['pore.test_bool'].size

    def test_stitch_radius_no_connections(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc * Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='radius', len_max=0,
                            label_stitches=['test', 'test2'])
        assert pn.Nt == (pn2.Nt * 2)

    def test_stitch_10_connections(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc * Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            label_stitches=['test', 'test2'])
        assert pn.Nt == (pn2.Nt * 2 + 10)

    def test_stitch_with_multiple_labels(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc * Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            label_stitches=['test', 'test2'])
        assert 'throat.test' in pn.keys()
        assert 'throat.test2' in pn.keys()

    def test_stitch_repeatedly(self):
        pn = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn2 = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn2['pore.coords'] += [1e-4 * 10, 0, 0]
        pn3 = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn3['pore.coords'] += [1e-4 * 20, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='nearest')
        op.topotools.stitch(network=pn, donor=pn3,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='nearest')
        assert pn.Nt == (pn2.Nt * 3 + 20)

    def test_dimensionality(self):
        # 3D network
        pn = op.network.Cubic(shape=[3, 4, 5])
        dims = op.topotools.dimensionality(pn)
        assert np.allclose(dims, np.array([True, True, True]))
        # 2D network
        pn = op.network.Cubic(shape=[3, 1, 5])
        dims = op.topotools.dimensionality(pn)
        assert np.allclose(dims, np.array([True, False, True]))
        # 1D network
        pn = op.network.Cubic(shape=[1, 1, 5], spacing=1.2345e-30)
        pn["pore.coords"][0, 0] = pn["pore.coords"][0, :].mean()
        dims = op.topotools.dimensionality(pn)
        assert np.allclose(dims, np.array([False, False, True]))

    def test_is_fully_connected(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        op.topotools.trim(network=pn, pores=[1, 5, 9, 13])
        assert op.topotools.is_fully_connected(pn) is False
        assert op.topotools.is_fully_connected(pn, pores_BC=[0, 1, 9, 10]) is True
        assert op.topotools.is_fully_connected(pn, pores_BC=[0, 10]) is True
        assert op.topotools.is_fully_connected(pn, pores_BC=[1, 10]) is False
        assert op.topotools.is_fully_connected(pn, pores_BC=[0, 9]) is False
        assert op.topotools.is_fully_connected(pn, pores_BC=[0]) is False
        assert op.topotools.is_fully_connected(pn, pores_BC=[]) is False

    def test_is_fully_connected_with_alg(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        Ts = pn.find_neighbor_throats(5)
        op.topotools.trim(network=pn, throats=Ts)
        assert op.topotools.is_fully_connected(pn) is False
        phase = op.phase.Phase(network=pn)
        phase['throat.diffusive_conductance'] = 1.0
        alg = op.algorithms.FickianDiffusion(network=pn, phase=phase)
        alg.set_value_BC(pores=[0, 1, 2, 3], values=1.0)
        alg.set_value_BC(pores=[12, 13, 14, 15], values=0.5)
        with pytest.raises(Exception):
            alg.run()
        alg.set_value_BC(pores=5, values=1.0)
        alg.run()
        assert np.ptp(alg['pore.concentration']) == 0.5

    def test_spacing_could_not_be_found(self):
        net = op.network.Cubic(shape=[1, 5, 1], spacing=1)
        net["pore.coords"][4, 1] += 5
        with pytest.raises(Exception):
            _ = net.spacing

    def test_spacing_1D(self):
        net = op.network.Cubic(shape=[2, 1, 1], spacing=1)
        assert np.all(op.topotools.get_spacing(net) == [1.0, 0.0, 0.0])

    def test_spacing_2D(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 1])
        assert np.all(op.topotools.get_spacing(net) == [1.0, 1.0, 0.0])

    def test_spacing_3D(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
        assert np.all(op.topotools.get_spacing(net) == [1.0, 1.0, 1.0])

    def test_spacing_2D_uneven(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 2])
        assert np.all(op.topotools.get_spacing(net) == [1.0, 2.0, 0.0])

    def test_spacing_3D_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        assert np.all(op.topotools.get_spacing(net) == [1.0, 2.0, 3.0])

    def test_shape_1D(self):
        net = op.network.Cubic(shape=[5, 1, 1])
        assert np.all(op.topotools.get_shape(net) == [5, 1, 1])
        net = op.network.Cubic(shape=[1, 1, 5])
        assert np.all(op.topotools.get_shape(net) == [1, 1, 5])
        net = op.network.Cubic(shape=[1, 5, 1])
        assert np.all(op.topotools.get_shape(net) == [1, 5, 1])
        net = op.network.Cubic(shape=5)
        assert np.all(op.topotools.get_shape(net) == [5, 1, 1])
        net = op.network.Cubic(shape=[5])
        assert np.all(op.topotools.get_shape(net) == [5, 1, 1])

    def test_shape_2D(self):
        net = op.network.Cubic(shape=[5, 5, 1])
        assert np.all(op.topotools.get_shape(net) == [5, 5, 1])
        net = op.network.Cubic(shape=[5, 1, 5])
        assert np.all(op.topotools.get_shape(net) == [5, 1, 5])
        net = op.network.Cubic(shape=[1, 5, 5])
        assert np.all(op.topotools.get_shape(net) == [1, 5, 5])
        net = op.network.Cubic(shape=[5, 2])
        assert np.all(op.topotools.get_shape(net) == [5, 2, 1])

    def test_shape_3D(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        assert np.all(op.topotools.get_shape(net) == [5, 5, 5])
        net = op.network.Cubic(shape=[3, 4, 5])
        assert np.all(op.topotools.get_shape(net) == [3, 4, 5])
        net = op.network.Cubic(shape=[1, 5, 1])
        assert np.all(op.topotools.get_shape(net) == [1, 5, 1])

    def test_shape_extremely_small_spacing(self):
        net = op.network.Cubic(shape=[200, 200, 1], spacing=1e-9)
        assert np.allclose(op.topotools.get_shape(net), [200, 200, 1])

    def test_spacing_3d_rotated(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[1, 1, 1])
        theta = 0.1
        R = np.array([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]], dtype=float)
        net['pore.coords'] = net['pore.coords'] @ R.T
        assert np.allclose(op.topotools.get_spacing(net), [1, 1, 1])

    def test_spacing_3d_rotated_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        theta = 0.1
        R = np.array([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])
        net['pore.coords'] = np.tensordot(net['pore.coords'], R, axes=(1, 1))
        assert np.allclose(op.topotools.get_spacing(net), [1, 2, 3])

    def test_spacing_2d_sheared(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=1)
        S = np.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert np.allclose(op.topotools.get_spacing(net), [1, 2**0.5, 0])

    def test_spacing_2d_sheared_uneven(self):
        net = op.network.Cubic(shape=[5, 5, 1], spacing=[1, 2])
        S = np.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert np.allclose(op.topotools.get_spacing(net), [1.0, 2*(2**0.5), 0.0])

    def test_spacing_3d_sheared(self):
        net = op.network.Cubic(shape=[5, 5, 3], spacing=1)
        S = np.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert np.allclose(op.topotools.get_spacing(net), [1.0, 2**0.5, 1.0])

    def test_spacing_3d_sheared_uneven(self):
        net = op.network.Cubic(shape=[3, 4, 5], spacing=[1, 2, 3])
        S = np.array([[1, 1, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        net['pore.coords'] = (S@net['pore.coords'].T).T
        assert np.allclose(op.topotools.get_spacing(net), [1, 2*(2**0.5), 3])

    def test_spacing_on_joggled_network(self):
        net = op.network.Cubic(shape=[3, 4, 5])
        net['pore.coords'] += np.random.rand(net.Np, 3)
        with pytest.raises(Exception):
            _ = op.topotools.get_spacing(net)

    def test_spacing_on_network_with_boundary_pores(self):
        net = op.network.Cubic(shape=[3, 4, 5])
        net.add_boundary_pores()
        with pytest.raises(Exception):
            _ = op.topotools.get_spacing(net)

    def test_find_interface_throats(self):
        net = op.network.Cubic([5, 1, 1])
        with pytest.raises(Exception):
            op.topotools.find_interface_throats(net, [1, 2, 3], [3, 4, 5])
        # Test pores that do have an interface
        Ts_int = op.topotools.find_interface_throats(net, 0, 1)
        assert_allclose(Ts_int, np.array([0]))
        # Test pores without an interface
        Ts_int = op.topotools.find_interface_throats(net, 2, 4)
        assert_allclose(Ts_int, np.array([]))
        # Test calling the function with labels
        net["pore.internal"] = ~net["pore.surface"]
        P1 = net.pores('internal')
        P2 = net.pores('surface')
        Ts_int = op.topotools.find_interface_throats(net, P1, P2)
        assert_allclose(Ts_int, np.array([0, 3]))

    def test_get_domain_area(self):
        net = op.network.Cubic([5, 4, 3])
        area = op.topotools.get_domain_area(
            network=net,
            inlets=net.pores('left'),
            outlets=net.pores('right')
        )
        assert area == 6.0
        area = op.topotools.get_domain_area(
            network=net,
            inlets=net.pores('top'),
            outlets=net.pores('bottom')
        )
        assert area == 12.0
        area = op.topotools.get_domain_area(
            network=net,
            inlets=net.pores('front'),
            outlets=net.pores('back')
        )
        assert area == 8.0

    def test_get_domain_length(self):
        net = op.network.Cubic([5, 4, 3])
        length = op.topotools.get_domain_length(
            network=net,
            inlets=net.pores('left'),
            outlets=net.pores('right')
        )
        assert length == 4.0
        length = op.topotools.get_domain_length(
            network=net,
            inlets=net.pores('top'),
            outlets=net.pores('bottom')
        )
        assert length == 2.0
        length = op.topotools.get_domain_length(
            network=net,
            inlets=net.pores('front'),
            outlets=net.pores('back')
        )
        assert length == 3.0


if __name__ == '__main__':

    t = TopotoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
