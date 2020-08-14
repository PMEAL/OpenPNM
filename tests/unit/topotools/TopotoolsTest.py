import pytest
import numpy as np
import openpnm as op
from numpy.testing import assert_allclose
from openpnm import topotools


class TopotoolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_reduce_coordination(self):
        net = op.network.Cubic(shape=[10, 10, 10], connectivity=26)
        a = np.mean(net.num_neighbors(pores=net.Ps, flatten=False))
        b = 20.952
        assert a == b
        topotools.reduce_coordination(network=net, z=6)
        a = np.mean(net.num_neighbors(pores=net.Ps, flatten=False))
        b = 6.0
        assert_allclose(a, b)
        h = net.check_network_health()
        assert h.health

    def test_label_faces(self):
        net = op.network.Cubic(shape=[3, 3, 3], connectivity=6)
        net.clear(mode='labels')
        assert net.labels() == ['pore.all', 'throat.all']
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
        assert net.labels() == ['pore.all', 'throat.all']
        topotools.find_surface_pores(network=net)
        assert net.num_pores('surface') == 66

    def test_find_surface_pores_custom_markers_2d(self):
        net = op.network.Cubic(shape=[4, 4, 1], spacing=1)
        net.clear(mode='labels')
        assert net.labels() == ['pore.all', 'throat.all']
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
        assert net.labels() == ['pore.all', 'throat.all']
        markers = [[-1, 2, 2], [2, -1, 2], [2, 5, 2], [5, 2, 2]]
        topotools.find_surface_pores(network=net, markers=markers)
        assert net.num_pores('surface') == 48
        markers = [[-1], [2], [2], [5]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)
        markers = [[-1, 2], [2, -1], [2, 5], [5, 2]]
        with pytest.raises(Exception):
            topotools.find_surface_pores(network=net, markers=markers)

    def test_find_pore_to_pore_distance(self):
        net = op.network.Cubic(shape=[3, 3, 3], connectivity=6)
        dm = topotools.find_pore_to_pore_distance(network=net,
                                                  pores1=net.pores('left'),
                                                  pores2=net.pores('right'))
        a = np.unique(dm)
        b = np.array([2., 2.23606798, 2.44948974, 2.82842712, 3., 3.46410162])
        assert np.allclose(a, b)

    def test_template_sphere_shell(self):
        im = topotools.template_sphere_shell(outer_radius=4, inner_radius=2)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 218
        assert net.Nt == 480
        im = topotools.template_sphere_shell(outer_radius=4)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 251
        assert net.Nt == 618

    def test_template_cylinder_annulus(self):
        im = topotools.template_cylinder_annulus(height=10, outer_radius=4,
                                                 inner_radius=2)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 320
        assert net.Nt == 688
        im = topotools.template_cylinder_annulus(height=10, outer_radius=4)
        net = op.network.CubicTemplate(template=im, spacing=1)
        assert net.Np == 450
        assert net.Nt == 1165
        im = topotools.template_cylinder_annulus(height=1, outer_radius=4)
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
        pn = op.network.Cubic(shape=[3, 3, 3])
        pn2 = op.network.Cubic(shape=[3, 3, 3])
        pn2['pore.coords'] += [0, 0, 3]
        pn2['pore.net_02'] = True
        pn2['throat.net_02'] = True
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
        geo2 = op.geometry.StickAndBall(network=pn2, pores=pn2.Ps, throats=pn2.Ts)
        geo2['pore.test_vals'] = 1.5
        op.topotools.merge_networks(network=pn, donor=pn2)
        assert isinstance(pn['pore.test_vals'], np.ndarray)
        assert ('pore.' + geo2.name) in pn.keys()

    def test_subdivide_3D(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        assert net.Np == 27
        assert net.Nt == 54
        op.topotools.subdivide(net, pores=13, shape=[5, 5, 5], labels="blah")
        assert net.pores("blah").size == 125
        assert net.throats("blah").size == 300 + 25 * 6
        assert net.Np == 27 - 1 + 125
        assert net.Nt == 54 - 6 + 300 + 25 * 6

    def test_subdivide_2D(self):
        net = op.network.Cubic(shape=[1, 3, 3])
        assert net.Np == 9
        assert net.Nt == 12
        op.topotools.subdivide(net, pores=4, shape=[1, 5, 5], labels="blah")
        assert net.pores("blah").size == 25
        assert net.throats("blah").size == 40 + 5 * 4
        assert net.Np == 9 - 1 + 25
        assert net.Nt == 12 - 4 + 40 + 5 * 4

    def test_merge_pores(self):
        testnet = op.network.Cubic(shape=[10, 10, 10])
        to_merge = [[0, 1], [998, 999]]
        topotools.merge_pores(testnet, to_merge)
        assert testnet.Np == 998

    def test_merge_pores_coords(self):
        r"""
        Coordinates of merged pores should be centroid of the enclosing convex
        hull.

        This test verifies that if one subdivides a pore and then merge it
        with a bunch of other pores, the coordinates of the new pore should be
        exactly the same as when one merges the same pores without subdiving.

        """
        # Subdivide first, then merge
        testnet = op.network.Cubic(shape=[1, 10, 10])
        testnet["pore.to_merge"] = False
        testnet["pore.to_merge"][[14, 15, 16, 24, 25, 26, 34, 35, 36]] = True
        topotools.subdivide(testnet, pores=15, shape=[1, 10, 10],
                            labels="subdivided")
        topotools.merge_pores(testnet, labels="new_pore",
                              pores=testnet.pores(["subdivided", "to_merge"]))
        xyz_w_subdivide = testnet['pore.coords'][testnet.pores("new_pore")]

        # No subdivide, only merge
        testnet = op.network.Cubic(shape=[1, 10, 10])
        testnet["pore.to_merge"] = False
        testnet["pore.to_merge"][[14, 15, 16, 24, 25, 26, 34, 35, 36]] = True
        topotools.merge_pores(testnet, labels="new_pore",
                              pores=testnet.pores("to_merge"))
        xyz_wo_subdivide = testnet['pore.coords'][testnet.pores("new_pore")]

        # Compare the two coords
        assert_allclose(xyz_w_subdivide, xyz_wo_subdivide)

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

    def test_ispercolating(self):
        net = op.network.Cubic(shape=[10, 10, 10], connectivity=26)
        tmask = net['throat.all']
        Pin = net.pores('left')
        Pout = net.pores('right')
        am = net.create_adjacency_matrix(weights=tmask, fmt='coo')
        val = topotools.ispercolating(am=am, mode='bond',
                                      inlets=Pin, outlets=Pout)
        assert val
        val = topotools.ispercolating(am=am, mode='site',
                                      inlets=Pin, outlets=Pout)
        assert val

    def test_trim_pores(self):
        np.random.seed(1)
        pn = op.network.Cubic(shape=[2, 2, 2], spacing=1)
        Ps = pn.pores()[:4]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        geo1 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
        Ps = pn.pores()[4:]
        Ts = pn.find_neighbor_throats(pores=Ps, mode='union')
        geo2 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
        geo1['pore.random'] = np.random.random(geo1.Np)
        geo2['pore.random'] = np.random.random(geo2.Np)
        trimmers = pn['pore.random'] < 0.25
        topotools.trim(pn, pores=pn.pores()[trimmers])
        assert ~np.any(pn['pore.random'] < 0.25)

    def test_trim_throats(self):
        np.random.seed(1)
        pn = op.network.Cubic(shape=[2, 2, 2], spacing=5)
        Ps = pn.pores()[:4]
        Ts1 = pn.find_neighbor_throats(pores=Ps, mode='or')
        geo1 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts1)
        Ps = pn.pores()[4:]
        Ts2 = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        geo2 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts2)
        geo1['throat.random'] = np.random.random(geo1.Nt)
        geo2['throat.random'] = np.random.random(geo2.Nt)
        trimmers = pn['throat.random'] < 0.25
        topotools.trim(pn, throats=pn.throats()[trimmers])
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
        op.topotools.extend(network=pn, pore_coords=[[3, 3, 3], [3, 3, 4]])
        assert np.any(np.isnan(pn['pore.test_float']))
        assert np.any(np.isnan(pn['pore.test_int']))
        assert pn['pore.test_bool'].sum() < pn['pore.test_bool'].size

    def test_extend_geometry_present(self):
        pn = op.network.Cubic(shape=[2, 2, 1])
        geo = op.geometry.StickAndBall(network=pn)
        geo['pore.test_float'] = 1.0
        geo['pore.test_int'] = 1
        geo['pore.test_bool'] = True
        op.topotools.extend(network=pn, pore_coords=[[3, 3, 3], [3, 3, 4]])
        assert ~np.any(np.isnan(geo['pore.test_float']))
        assert ~np.any(np.isnan(geo['pore.test_int']))
        assert geo['pore.test_bool'].sum() == geo['pore.test_bool'].size

    def test_extend_phase_present(self):
        pn = op.network.Cubic(shape=[2, 2, 1])
        air = op.phases.Air(network=pn)
        air['pore.test_float'] = 1.0
        air['pore.test_int'] = 1
        air['pore.test_bool'] = True
        with pytest.raises(Exception):
            op.topotools.extend(network=pn, pore_coords=[[3, 3, 3], [3, 3, 4]])

    def test_stitch_radius_no_connections(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc*Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='radius', len_max=0,
                            label_stitches=['test', 'test2'])
        assert pn.Nt == (pn2.Nt*2)

    def test_stitch_10_connections(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc*Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            label_stitches=['test', 'test2'])
        assert pn.Nt == (pn2.Nt*2 + 10)

    def test_stitch_with_multiple_labels(self):
        Nx, Ny, Nz = (10, 10, 1)
        Lc = 1e-4
        pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2 = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
        pn2['pore.coords'] += [Lc*Nx, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            label_stitches=['test', 'test2'])
        assert 'throat.test' in pn.keys()
        assert 'throat.test2' in pn.keys()

    def test_stitch_repeatedly(self):
        pn = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn2 = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn2['pore.coords'] += [1e-4*10, 0, 0]
        pn3 = op.network.Cubic(shape=[10, 10, 1], spacing=1e-4)
        pn3['pore.coords'] += [1e-4*20, 0, 0]
        op.topotools.stitch(network=pn, donor=pn2,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='nearest')
        op.topotools.stitch(network=pn, donor=pn3,
                            P_network=pn.pores('back'), P_donor=pn2.pores('front'),
                            method='nearest')
        assert pn.Nt == (pn2.Nt*3 + 20)

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
        pn = op.network.Cubic(shape=[1, 1, 5])
        dims = op.topotools.dimensionality(pn)
        assert np.allclose(dims, np.array([False, False, True]))


if __name__ == '__main__':

    t = TopotoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: ' + item)
            t.__getattribute__(item)()
