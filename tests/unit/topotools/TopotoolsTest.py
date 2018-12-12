import openpnm as op
import numpy as np
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

    def test_find_surface_pores(self):
        from skimage.morphology import ball
        net = op.network.CubicTemplate(template=ball(3), spacing=1)
        net.clear(mode='labels')
        assert net.labels() == ['pore.all', 'throat.all']
        topotools.find_surface_pores(network=net)
        assert net.num_pores('surface') == 42

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
        net2['pore.test3'] = True
        net2['pore.test4'] = 10.0
        topotools.merge_networks(net1, net2)
        assert np.sum(net1['pore.test1']) == 27
        assert np.sum(net1['pore.test3']) == 27
        assert np.sum(net1['pore.test2'][:27]) == 270
        assert np.sum(net1['pore.test4'][27:]) == 270
        assert 'pore.test1' not in net2
        assert 'pore.test2' not in net2

    def test_merge_pores(self):
        testnet = op.network.Cubic(shape=[10, 10, 10])
        xyz_old = testnet['pore.coords'].copy()
        to_merge = [[0, 1], [998, 999]]
        topotools.merge_pores(testnet, to_merge)
        xyz = testnet['pore.coords']
        xyz1 = xyz[-2]
        xyz2 = xyz[-1]
        xyz1_desired = xyz_old[0:2].mean(axis=0)
        xyz2_desired = xyz_old[998::].mean(axis=0)
        assert testnet.Np == 998
        assert_allclose(xyz1, xyz1_desired)
        assert_allclose(xyz2, xyz2_desired)

    def test_connect_pores(self):
        testnet = op.network.Cubic(shape=[10, 10, 10])
        Nt_old= testnet.Nt
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


if __name__ == '__main__':

    t = TopotoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
