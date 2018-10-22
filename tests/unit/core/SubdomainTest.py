import openpnm as op
import scipy as sp
import pytest


class SubdomainTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = sp.rand(self.net.Np)
        self.geo.add_model(propname='pore.volume',
                           model=op.models.geometry.pore_volume.sphere)
        self.geo['throat.diameter'] = sp.rand(self.net.Nt)
        self.geo.add_model(propname='throat.area',
                           model=op.models.geometry.throat_area.cylinder)
        self.geo.regenerate_models()
        self.phase1 = op.phases.GenericPhase(network=self.net)
        self.phase2 = op.phases.GenericPhase(network=self.net)
        self.phys1 = op.physics.GenericPhysics(network=self.net,
                                               geometry=self.geo,
                                               phase=self.phase1)
        self.phys1['pore.blah'] = 1.0
        self.phys2 = op.physics.GenericPhysics(network=self.net,
                                               geometry=self.geo,
                                               phase=self.phase2)
        self.phys2['pore.blah'] = 2.0

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_drop_locations_from_geom_successively_with_single_geometry(self):
        assert self.geo.Np == 27
        assert self.geo.Nt == 54
        self.geo.drop_locations(pores=[0, 1, 2], throats=[0, 1, 2])
        assert self.geo.Np == 24
        assert self.geo.Nt == 51
        self.geo.drop_locations(pores=[3, 4], throats=[3, 4])
        assert self.geo.Np == 22
        assert self.geo.Nt == 49
        self.geo.add_locations(pores=[0, 1, 2, 3, 4], throats=[0, 1, 2, 3, 4])
        assert self.geo.Np == 27
        assert self.geo.Nt == 54

    def test_drop_locations_from_physics_successively_with_two_physics(self):
        assert self.phys1.Np == 27
        assert self.phys1.Nt == 54
        self.phys1.drop_locations(pores=[0, 1], throats=[0, 1])
        assert self.phys1.Np == 25
        assert self.phys1.Nt == 52
        self.phys1.drop_locations(pores=[3, 4], throats=[3, 4])
        assert self.phys1.Np == 23
        assert self.phys1.Nt == 50
        self.phys1.add_locations(pores=[0, 1, 3, 4], throats=[0, 1, 3, 4])
        assert self.phys1.Np == 27
        assert self.phys1.Nt == 54

    def test_drop_locations_all_but_not_complete(self):
        assert self.phys1.Np == 27
        assert self.phys1.Nt == 54
        assert 'pore.'+self.phys1.name in self.phase1.keys()
        assert 'throat.'+self.phys1.name in self.phase1.keys()
        self.phys1.drop_locations(pores=self.net.Ps)
        assert 'pore.'+self.phys1.name in self.phase1.keys()
        assert self.phase1.num_pores(self.phys1.name) == 0
        assert 'throat.'+self.phys1.name in self.phase1.keys()
        self.phys1.drop_locations(throats=self.net.Ts)
        assert 'throat.'+self.phys1.name in self.phase1.keys()
        assert self.phase1.num_throats(self.phys1.name) == 0
        self.phys1.add_locations(pores=self.net.Ps, throats=self.net.Ts)

    def test_writting_subdict_names_across_subdomains(self):
        ws = op.Workspace()
        proj = ws.new_project()

        pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
        Ps = pn['pore.coords'][:, 0] < pn['pore.coords'][:, 0].mean()
        Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        geo1 = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)

        Ps = pn['pore.coords'][:, 0] >= pn['pore.coords'][:, 0].mean()
        Ts = pn.find_neighbor_throats(pores=Ps, mode='or')
        geo2 = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)

        pn['pore.foo'] = 1
        # Can't create a subdict below foo
        with pytest.raises(Exception):
            pn['pore.foo.bar'] = 1
        # Can create a subdict directly
        pn['pore.baz.bar'] = 2
        # Can't create a new item already used as subdict
        with pytest.raises(Exception):
            pn['pore.baz'] = 2

        # Also works on subdomains
        geo1['pore.blah'] = 1
        with pytest.raises(Exception):
            geo1['pore.blah.boo'] = 1
        geo1['pore.bee.bop'] = 1
        with pytest.raises(Exception):
            geo1['pore.bee'] = 1

        # Now start looking across objects
        with pytest.raises(Exception):
            geo1['pore.foo'] = 1  # Already exists on pn
        with pytest.raises(Exception):
            geo1['pore.foo.bar'] = 1  # pore.foo already exists on pn
        with pytest.raises(Exception):
            geo1['pore.baz'] = 1  # pore.baz.bar already exists on pn

        # Now start looking across objects
        geo2['pore.blah'] = 1
        geo2['pore.bee.bop'] = 1
        with pytest.raises(Exception):
            geo1['pore.bee'] = 1

        with pytest.raises(Exception):
            pn['pore.bee'] = 1

        with pytest.raises(Exception):
            pn['pore.bee.bop'] = 1


if __name__ == '__main__':

    t = SubdomainTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
