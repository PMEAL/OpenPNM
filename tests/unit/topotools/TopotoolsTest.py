import openpnm as op
import numpy as np
import pytest
from numpy.testing import assert_approx_equal
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
        assert_approx_equal(a, b)
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
        with pytest.warns(UserWarning):
            topotools.merge(net1, net2)
        assert np.sum(net1['pore.test1']) == 27
        assert np.sum(net1['pore.test3']) == 27
        assert np.sum(net1['pore.test2']) == 270
        assert np.sum(np.isnan(net1['pore.test4'])) == 27


if __name__ == '__main__':

    t = TopotoolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
