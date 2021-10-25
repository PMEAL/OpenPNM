import pytest
import numpy as np
import openpnm as op


class SubdomainTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = np.random.rand(self.net.Np)
        self.geo.add_model(propname='pore.volume',
                           model=op.models.geometry.pore_volume.sphere)
        self.geo['throat.diameter'] = np.random.rand(self.net.Nt)
        self.geo.add_model(propname='throat.area',
                           model=op.models.geometry.throat_cross_sectional_area.cylinder)
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

    def test_interleaving_missing_objects(self):
        pn = op.network.Cubic(shape=[3, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=pn, pores=[0], throats=[0])
        geo2 = op.geometry.GenericGeometry(network=pn, pores=[1], throats=[1])

        geo1['pore.test_float'] = 1.0
        geo2['pore.test_float'] = 2.0
        assert pn['pore.test_float'].dtype == float
        assert np.isnan(pn['pore.test_float']).sum() == 1

        geo1['pore.test_int'] = 1
        geo2['pore.test_int'] = 2
        assert pn['pore.test_int'].dtype == float
        assert np.isnan(pn['pore.test_int']).sum() == 1

        geo1['pore.test_bool'] = True
        geo2['pore.test_bool'] = False
        assert pn['pore.test_bool'].dtype == bool
        assert pn['pore.test_bool'].sum() == 1

        geo1['pore.test_int'] = 1.0
        geo2['pore.test_int'] = 2
        assert pn['pore.test_int'].dtype == float
        assert np.isnan(pn['pore.test_int']).sum() == 1

        geo1['pore.test_int'] = 1
        geo2['pore.test_int'] = 2.0
        assert pn['pore.test_int'].dtype == float
        assert np.isnan(pn['pore.test_int']).sum() == 1

    def test_interleaving_mixed_data(self):
        pn = op.network.Cubic(shape=[3, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=pn, pores=[0], throats=[0])
        geo2 = op.geometry.GenericGeometry(network=pn, pores=[1], throats=[1])
        geo3 = op.geometry.GenericGeometry(network=pn, pores=[2])

        # Test floats with all arrays present
        geo1['pore.test_float'] = 1.0
        geo2['pore.test_float'] = 2.0
        geo3['pore.test_float'] = 3.0
        assert pn['pore.test_float'].dtype == float
        assert np.isnan(pn['pore.test_float']).sum() == 0

        # Test mixed datatype with all arrays present
        # It's not clear that we want this behavior
        # geo1['pore.test_mixed'] = 1.0
        # geo2['pore.test_mixed'] = 2
        # geo3['pore.test_mixed'] = False
        # assert pn['pore.test_mixed'].dtype == float
        # assert np.isnan(pn['pore.test_mixed']).sum() == 0

        # Check heterogeneous datatypes
        geo1['pore.test_mixed'] = 1
        geo2['pore.test_mixed'] = 2
        geo3['pore.test_mixed'] = 3.0
        assert pn['pore.test_mixed'].dtype == float

        # Make sure order doesn't matter
        geo1['pore.test_mixed'] = 1.0
        geo2['pore.test_mixed'] = 2
        geo3['pore.test_mixed'] = 3
        assert pn['pore.test_mixed'].dtype == float

    def test_interleaving_partial_data(self):
        pn = op.network.Cubic(shape=[3, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=pn, pores=[0], throats=[0])
        geo2 = op.geometry.GenericGeometry(network=pn, pores=[1], throats=[1])
        geo3 = op.geometry.GenericGeometry(network=pn, pores=[2])

        # Test ints with a missing array
        geo1['pore.test_int_missing'] = 1
        geo2['pore.test_int_missing'] = 2
        assert np.isnan(pn['pore.test_int_missing']).sum() == 1
        assert pn['pore.test_int_missing'].dtype == float

        # Test ints with all arrays present
        geo1['pore.test_int'] = 1
        geo2['pore.test_int'] = 2
        geo3['pore.test_int'] = 3
        assert pn['pore.test_int'].dtype == int

        # Test booleans with a missing array
        geo1['pore.test_bool'] = True
        geo2['pore.test_bool'] = False
        assert pn['pore.test_bool'].dtype == bool
        assert pn['pore.test_bool'].sum() == 1

    def test_interleave_data_bool(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        # Ensure Falses return in missing places
        geom1['pore.blah'] = True
        assert np.all(~geom2['pore.blah'])
        assert np.sum(net['pore.blah']) == 4
        # Ensure all Trues returned now
        geom2['pore.blah'] = True
        assert np.all(geom2['pore.blah'])
        assert np.sum(net['pore.blah']) == 8

    def test_interleave_data_int(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = 1
        # Ensure ints are returned geom1
        assert 'int' in geom1['pore.blah'].dtype.name
        # Ensure nans are returned on geom2
        assert np.all(np.isnan(geom2['pore.blah']))
        # Ensure interleaved array is float with nans
        assert 'float' in net['pore.blah'].dtype.name
        # Ensure missing values are floats
        assert np.sum(np.isnan(net['pore.blah'])) == 4

    def test_interleave_data_float(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = 1.0
        # Ensure flaots are returned geom1
        assert 'float' in geom1['pore.blah'].dtype.name
        # Ensure nans are returned on geom2
        assert np.all(np.isnan(geom2['pore.blah']))
        # Ensure interleaved array is float with nans
        assert 'float' in net['pore.blah'].dtype.name
        # Ensure missing values are floats
        assert np.sum(np.isnan(net['pore.blah'])) == 4

    def test_interleave_data_object(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        _ = op.geometry.GenericGeometry(network=net, pores=Ps)
        geom1['pore.blah'] = [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1]]
        assert 'object' in net['pore.blah'].dtype.name
        # Ensure missing elements are None
        assert np.sum([item is None for item in net['pore.blah']]) == 4

    def test_interleave_data_key_error(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        Ps = net.pores('top')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps)
        Ps = net.pores('bottom')
        geom2 = op.geometry.GenericGeometry(network=net, pores=Ps)
        with pytest.raises(KeyError):
            _ = net['pore.blah']
        with pytest.raises(KeyError):
            _ = geom1['pore.blah']
        with pytest.raises(KeyError):
            _ = geom2['pore.blah']

    def test_interleave_data_float_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = 1.0
        assert np.any(np.isnan(net['pore.blah']))

    def test_interleave_data_int_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = 1
        assert np.any(np.isnan(net['pore.blah']))

    def test_interleave_data_bool_missing_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geom = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2])
        geom['pore.blah'] = True
        assert np.sum(net['pore.blah']) == geom.Np

    def test_interleave_data_float_missing_physics(self):
        net = op.network.Cubic(shape=[4, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=net, pores=[0, 1], throats=[0, 1])
        op.geometry.GenericGeometry(network=net, pores=[2, 3], throats=[2])
        air = op.phases.Air(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo1)
        phys['pore.blah'] = 1.0
        assert np.any(np.isnan(air['pore.blah']))

    def test_interleave_data_int_missing_physics(self):
        net = op.network.Cubic(shape=[4, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=net, pores=[0, 1], throats=[0, 1])
        op.geometry.GenericGeometry(network=net, pores=[2, 3], throats=[2])
        air = op.phases.Air(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo1)
        phys['pore.blah'] = 1
        assert np.any(np.isnan(air['pore.blah']))

    def test_interleave_data_bool_missing_physics(self):
        net = op.network.Cubic(shape=[4, 1, 1])
        geo1 = op.geometry.GenericGeometry(network=net, pores=[0, 1], throats=[0, 1])
        op.geometry.GenericGeometry(network=net, pores=[2, 3], throats=[2])
        air = op.phases.Air(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=air, geometry=geo1)
        phys['pore.blah'] = True
        assert np.sum(air['pore.blah']) == phys.Np

    def test_writing_subdict_names_across_subdomains(self):
        ws = op.Workspace()
        proj = ws.new_project()

        pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
        Ps = pn['pore.coords'][:, 0] < pn['pore.coords'][:, 0].mean()
        Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
        geo1 = op.geometry._StickAndBall(network=pn, pores=Ps, throats=Ts)

        Ps = pn['pore.coords'][:, 0] >= pn['pore.coords'][:, 0].mean()
        Ts = pn.find_neighbor_throats(pores=Ps, mode='or')
        geo2 = op.geometry._StickAndBall(network=pn, pores=Ps, throats=Ts)

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
            print(f'Running test: {item}')
            t.__getattribute__(item)()
