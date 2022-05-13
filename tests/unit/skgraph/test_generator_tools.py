import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from openpnm._skgraph import generators as gen
from openpnm._skgraph import tools


class SKGRGeneratorToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_parse_points(self):
        pts = gen.tools.parse_points(shape=[1, 1, 1], points=10)
        assert pts.shape == (10, 3)
        assert pts[:, 2].max() > 0
        pts = gen.tools.parse_points(shape=[1, 1, 0], points=10)
        assert pts.shape == (10, 3)
        assert pts[:, 2].max() == 0
        pts = gen.tools.parse_points(shape=[1, 1], points=10)
        assert pts.shape == (10, 3)
        assert pts[:, 2].max() > 0
        r, q, z = tools.cart2cyl(*pts.T)
        assert r.max() < 1
        assert r.max() > 0
        pts = gen.tools.parse_points(shape=[1, 0], points=10)
        assert pts.shape == (10, 3)
        assert pts[:, 2].max() == 0
        r, q, z = tools.cart2cyl(*pts.T)
        assert r.max() < 1
        assert r.max() > 0
        pts = gen.tools.parse_points(shape=[1], points=10)
        assert pts.shape == (10, 3)
        assert pts[:, 2].max() > 0
        r, q, p = tools.cart2sph(*pts.T)
        assert r.max() < 1
        assert r.max() > 0

    def test_add_all_label(self):
        d = gen.cubic(shape=[3, 3, 3], node_prefix='pore', edge_prefix='throat')
        assert len(d.keys()) == 2
        d = gen.tools.add_all_label(d)
        assert len(d.keys()) == 4
        assert 'pore.all' in d.keys()
        assert 'throat.all' in d.keys()

    def test_label_faces_cubic(self):
        d = gen.cubic(shape=[3, 3, 3], node_prefix='pore', edge_prefix='throat')
        d = gen.tools.label_faces_cubic(d, rtol=0.1)
        assert d['pore.left'].sum() == 9
        assert d['pore.right'].sum() == 9
        assert d['pore.front'].sum() == 9
        assert d['pore.back'].sum() == 9
        assert d['pore.top'].sum() == 9
        assert d['pore.bottom'].sum() == 9

    def test_template_sphere_shell(self):
        im = gen.tools.template_sphere_shell(r_outer=10, r_inner=0)
        assert im.sum() == 4139
        im = gen.tools.template_sphere_shell(r_outer=10, r_inner=5)
        assert im.sum() == 3624

    def test_template_cylinder_annulus(self):
        im = gen.tools.template_cylinder_annulus(z=10, r_outer=10, r_inner=0)
        assert im.sum() == 3050
        im = gen.tools.template_cylinder_annulus(z=10, r_outer=10, r_inner=5)
        assert im.sum() == 2240
        im = gen.tools.template_cylinder_annulus(z=0, r_outer=10, r_inner=0)
        assert im.sum() == 305
        im = gen.tools.template_cylinder_annulus(z=0, r_outer=10, r_inner=5)
        assert im.sum() == 224

    def test_generate_base_points_and_reflect(self):
        f = gen.tools.generate_base_points
        pts = f(20, domain_size=[1, 1, 1], reflect=False)
        assert pts.shape == (20, 3)
        pts = f(20, domain_size=[1, 1, 1], reflect=True)
        assert pts.shape == (140, 3)
        pts = f(20, domain_size=[1, 1, 0], reflect=False)
        assert pts.shape == (20, 3)
        pts = f(20, domain_size=[1, 1, 0], reflect=True)
        assert pts.shape == (100, 3)
        pts = f(20, domain_size=[1, 1], reflect=False)
        assert pts.shape == (20, 3)
        pts = f(20, domain_size=[1, 1], reflect=True)
        assert pts.shape == (120, 3)
        pts = f(20, domain_size=[1, 0], reflect=False)
        assert pts.shape == (20, 3)
        pts = f(20, domain_size=[1, 0], reflect=True)
        assert pts.shape == (40, 3)
        pts = f(20, domain_size=[1], reflect=False)
        assert pts.shape == (20, 3)
        pts = f(20, domain_size=[1], reflect=True)
        assert pts.shape == (40, 3)


if __name__ == '__main__':

    t = SKGRGeneratorToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
