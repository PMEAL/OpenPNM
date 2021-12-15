import openpnm as op
from openpnm.io import Dict
import py
import os


class DictTest:

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
        Ts = self.net.find_neighbor_throats(pores=Ps, mode='xnor')
        self.geo_2 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_2['pore.boo'] = 1
        self.geo_2['throat.boo'] = 1

        self.phase_1 = op.phase.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phase.GenericPhase(network=self.net)
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
        ws = op.Workspace()
        ws.clear()

    def test_to_dict_missing_all_physics(self):
        net = op.network.Cubic(shape=[4, 4, 4])
        op.geometry.GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)
        phase = op.phase.GenericPhase(network=net)

        Dict.to_dict(network=net, phases=[phase], flatten=True,
                     interleave=True, categorize_by=[])

    def test_to_dict_flattened_interleaved(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=True, categorize_by=[])

        a = set(['net_01', 'phase_01', 'phase_02'])

        assert a == set(D.keys())
        assert set(['geo_01', 'geo_02']).isdisjoint(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase_02'].keys())

    def test_to_dict_flattened_not_interleaved(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False, categorize_by=[])

        a = set([i.name for i in self.net.project])

        assert a == set(D.keys())
        assert set(['geo_01', 'geo_02']).isdisjoint(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase_02'].keys())

    def test_to_dict_not_flattened_interleaved(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True, categorize_by=[])

        a = set(['net_01', 'phase_01', 'phase_02'])

        assert a == set(D.keys())

        assert set(['geo_01', 'geo_02']).isdisjoint(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase_02'].keys())

    def test_to_dict_not_flattened_not_interleaved(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False, categorize_by=[])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        _ = set(['pore', 'throat'])

        # Ensure NOT categorized by object
        assert b == set(D.keys())

        # Ensure NOT flattened
        assert set(['geo_01', 'geo_02']).issubset(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).issubset(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).issubset(D['phase_02'].keys())
        # Ensure no cross talk between phases
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase_02'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase_01'].keys())

    def test_to_dict_not_flat_not_interleaved_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['object'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        _ = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        _ = set(['pore', 'throat'])
        e = set(['network', 'phase'])

        # Ensure categorized by object
        assert e == set(D.keys())

        # Ensure flatten, which occurs when categorized by object
        keys = D['network']['net_01']['geometry'].keys()
        assert set(['geo_01', 'geo_02']).issubset(keys)
        keys = D['phase'].keys()
        assert set(['phase_01', 'phase_02']).issubset(keys)
        keys = D['phase']['phase_01']['physics'].keys()
        assert set(['phys_01', 'phys_02']).issubset(keys)

    def test_to_dict_not_flat_not_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['data'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        _ = set(['pore', 'throat'])

        # Ensure NOT categorized by object
        assert b == set(D.keys())

        # Ensure NOT flattened
        assert set(['geo_01', 'geo_02']).issubset(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).issubset(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).issubset(D['phase_02'].keys())

        # Ensure categorized by data
        assert c.issubset(D['net_01'].keys())
        assert c.issubset(D['phase_01'].keys())
        assert c.issubset(D['phase_02'].keys())
        assert c.issubset(D['net_01']['geo_01'].keys())
        assert c.issubset(D['net_01']['geo_02'].keys())
        assert c.issubset(D['phase_01']['phys_01'].keys())
        assert c.issubset(D['phase_01']['phys_02'].keys())
        assert c.issubset(D['phase_02']['phys_03'].keys())
        assert c.issubset(D['phase_02']['phys_04'].keys())

    def test_to_dict_not_flat_not_interleaved_categorized_by_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['element'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

        # Ensure NOT categorized by object
        assert b == set(D.keys())

        # Ensure NOT flattened
        assert set(['geo_01', 'geo_02']).issubset(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).issubset(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).issubset(D['phase_02'].keys())

        # Ensure it's categorized by element
        assert d.issubset(D['net_01'].keys())
        assert d.issubset(D['phase_01'].keys())
        assert d.issubset(D['phase_01'].keys())
        assert d.issubset(D['phase_02'].keys())
        assert d.issubset(D['phase_02'].keys())
        assert d.issubset(D['net_01']['geo_01'].keys())
        assert d.issubset(D['net_01']['geo_01'].keys())
        assert d.issubset(D['net_01']['geo_02'].keys())
        assert d.issubset(D['net_01']['geo_02'].keys())
        assert d.issubset(D['phase_01']['phys_01'].keys())
        assert d.issubset(D['phase_01']['phys_01'].keys())
        assert d.issubset(D['phase_01']['phys_02'].keys())
        assert d.issubset(D['phase_01']['phys_02'].keys())
        assert d.issubset(D['phase_02']['phys_03'].keys())
        assert d.issubset(D['phase_02']['phys_03'].keys())
        assert d.issubset(D['phase_02']['phys_04'].keys())
        assert d.issubset(D['phase_02']['phys_04'].keys())

    def test_to_dict_not_flat_not_interleaved_cat_by_element_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['element', 'data'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

        # Ensure NOT categorized by object
        assert b == set(D.keys())

        # Ensure NOT flattened
        assert set(['geo_01', 'geo_02']).issubset(D['net_01'].keys())
        assert set(['phys_01', 'phys_02']).issubset(D['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).issubset(D['phase_02'].keys())

        # Ensure categorized by data and element
        assert d.issubset(D['net_01']['properties'].keys())
        assert d.issubset(D['net_01']['labels'].keys())
        assert d.issubset(D['phase_01']['properties'].keys())
        assert d.issubset(D['phase_01']['labels'].keys())
        assert d.issubset(D['phase_02']['properties'].keys())
        assert d.issubset(D['phase_02']['labels'].keys())
        assert d.issubset(D['net_01']['geo_01']['properties'].keys())
        assert d.issubset(D['net_01']['geo_01']['labels'].keys())
        assert d.issubset(D['net_01']['geo_02']['properties'].keys())
        assert d.issubset(D['net_01']['geo_02']['labels'].keys())
        assert d.issubset(D['phase_01']['phys_01']['properties'].keys())
        assert d.issubset(D['phase_01']['phys_01']['labels'].keys())
        assert d.issubset(D['phase_01']['phys_02']['properties'].keys())
        assert d.issubset(D['phase_01']['phys_02']['labels'].keys())
        assert d.issubset(D['phase_02']['phys_03']['properties'].keys())
        assert d.issubset(D['phase_02']['phys_03']['labels'].keys())
        assert d.issubset(D['phase_02']['phys_04']['properties'].keys())
        assert d.issubset(D['phase_02']['phys_04']['labels'].keys())

    def test_to_dict_not_flat_not_interleaved_cat_by_element_data_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['element', 'data', 'object'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        _ = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        d = set(['pore', 'throat'])
        e = set(['network', 'phase'])

        # Check if categorized by object, but not flattened
        assert e == set(D.keys())
        assert 'geometry' in D['network']['net_01'].keys()
        assert 'physics' in D['phase']['phase_01'].keys()
        assert 'physics' in D['phase']['phase_02'].keys()

        # Ensure it's categorized by object, data, and element
        assert d.issubset(D['network']['net_01']['labels'].keys())
        assert d.issubset(D['phase']['phase_01']['properties'].keys())
        assert d.issubset(D['phase']['phase_01']['labels'].keys())
        assert d.issubset(D['phase']['phase_02']['properties'].keys())
        assert d.issubset(D['phase']['phase_02']['labels'].keys())
        path = D['network']['net_01']['geometry']['geo_01']['properties']
        assert d.issubset(path.keys())
        path = D['network']['net_01']['geometry']['geo_01']['labels']
        assert d.issubset(path.keys())
        path = D['network']['net_01']['geometry']['geo_02']['properties']
        assert d.issubset(path.keys())
        path = D['network']['net_01']['geometry']['geo_02']['labels']
        assert d.issubset(path.keys())
        path = D['phase']['phase_01']['physics']['phys_01']['properties']
        assert d.issubset(path.keys())
        path = D['phase']['phase_01']['physics']['phys_01']['labels']
        assert d.issubset(path.keys())
        path = D['phase']['phase_01']['physics']['phys_02']['properties']
        assert d.issubset(path.keys())
        path = D['phase']['phase_01']['physics']['phys_02']['labels']
        assert d.issubset(path.keys())
        path = D['phase']['phase_02']['physics']['phys_03']['properties']
        assert d.issubset(path.keys())
        path = D['phase']['phase_02']['physics']['phys_03']['labels']
        assert d.issubset(path.keys())
        path = D['phase']['phase_02']['physics']['phys_04']['properties']
        assert d.issubset(path.keys())
        path = D['phase']['phase_02']['physics']['phys_04']['labels']
        assert d.issubset(path.keys())

    def test_to_dict_not_flat_not_interleaved_cat_by_element_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['element', 'object'])

        _ = set(['network', 'phase', 'physics', 'geometry'])
        _ = set(['net_01', 'phase_01', 'phase_02'])
        _ = set(['labels', 'properties'])
        d = set(['pore', 'throat'])
        e = set(['network', 'phase'])

        # Check if categorized by object
        assert e == set(D.keys())

        # Check if categorized by object, but not flattened
        assert e == set(D.keys())
        assert 'geometry' in D['network']['net_01'].keys()
        assert 'physics' in D['phase']['phase_01'].keys()
        assert 'physics' in D['phase']['phase_02'].keys()

        # Ensure it's categorized by element
        assert d.issubset(D['network']['net_01'].keys())
        assert d.issubset(D['phase']['phase_01'].keys())
        assert d.issubset(D['phase']['phase_01'].keys())
        assert d.issubset(D['phase']['phase_02'].keys())
        assert d.issubset(D['phase']['phase_02'].keys())
        assert d.issubset(D['network']['net_01']['geometry']['geo_01'].keys())
        assert d.issubset(D['network']['net_01']['geometry']['geo_01'].keys())
        assert d.issubset(D['network']['net_01']['geometry']['geo_02'].keys())
        assert d.issubset(D['network']['net_01']['geometry']['geo_02'].keys())
        assert d.issubset(D['phase']['phase_01']['physics']['phys_01'].keys())
        assert d.issubset(D['phase']['phase_01']['physics']['phys_01'].keys())
        assert d.issubset(D['phase']['phase_01']['physics']['phys_02'].keys())
        assert d.issubset(D['phase']['phase_01']['physics']['phys_02'].keys())
        assert d.issubset(D['phase']['phase_02']['physics']['phys_03'].keys())
        assert d.issubset(D['phase']['phase_02']['physics']['phys_03'].keys())
        assert d.issubset(D['phase']['phase_02']['physics']['phys_04'].keys())
        assert d.issubset(D['phase']['phase_02']['physics']['phys_04'].keys())

    def test_to_dict_flat_not_interleaved_categorized_by_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['element'])

        assert set(D.keys()) == set([i.name for i in self.net.project])

        d = set(['pore', 'throat'])
        assert d.issubset(D['net_01'].keys())
        assert d.issubset(D['geo_01'].keys())
        assert d.issubset(D['geo_02'].keys())
        assert d.issubset(D['phase_01'].keys())
        assert d.issubset(D['phase_02'].keys())
        assert d.issubset(D['phys_01'].keys())
        assert d.issubset(D['phys_02'].keys())
        assert d.issubset(D['phys_03'].keys())
        assert d.issubset(D['phys_04'].keys())

    def test_to_dict_flat_not_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['data'])

        assert set(D.keys()) == set([i.name for i in self.net.project])

        c = set(['labels', 'properties'])
        assert c.issubset(D['net_01'].keys())
        assert c.issubset(D['geo_01'].keys())
        assert c.issubset(D['geo_02'].keys())
        assert c.issubset(D['phase_01'].keys())
        assert c.issubset(D['phase_02'].keys())
        assert c.issubset(D['phys_01'].keys())
        assert c.issubset(D['phys_02'].keys())
        assert c.issubset(D['phys_03'].keys())
        assert c.issubset(D['phys_04'].keys())

    def test_to_dict_flat_not_interleaved_categorized_by_data_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['data', 'element'])

        assert set(D.keys()) == set([i.name for i in self.net.project])

        d = set(['pore', 'throat'])
        assert d.issubset(D['net_01']['labels'].keys())
        assert d.issubset(D['net_01']['properties'].keys())
        assert d.issubset(D['geo_01']['labels'].keys())
        assert d.issubset(D['geo_01']['properties'].keys())
        assert d.issubset(D['geo_02']['labels'].keys())
        assert d.issubset(D['geo_02']['properties'].keys())
        assert d.issubset(D['phase_01']['labels'].keys())
        assert d.issubset(D['phase_01']['properties'].keys())
        assert d.issubset(D['phase_02']['labels'].keys())
        assert d.issubset(D['phase_02']['properties'].keys())
        assert d.issubset(D['phys_01']['labels'].keys())
        assert d.issubset(D['phys_01']['properties'].keys())
        assert d.issubset(D['phys_02']['labels'].keys())
        assert d.issubset(D['phys_02']['properties'].keys())
        assert d.issubset(D['phys_03']['labels'].keys())
        assert d.issubset(D['phys_03']['properties'].keys())
        assert d.issubset(D['phys_04']['labels'].keys())
        assert d.issubset(D['phys_04']['properties'].keys())

    def test_to_dict_interleaved_categorized_by_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['element'])

        b = set(['net_01', 'phase_01', 'phase_02'])
        assert set(D.keys()) == b

        d = set(['pore', 'throat'])
        assert d.issubset(D['net_01'].keys())
        assert d.issubset(D['phase_01'].keys())
        assert d.issubset(D['phase_02'].keys())

    def test_to_dict_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['data'])

        b = set(['net_01', 'phase_01', 'phase_02'])
        assert set(D.keys()) == b

        d = set(['labels', 'properties'])
        assert d.issubset(D['net_01'].keys())
        assert d.issubset(D['phase_01'].keys())
        assert d.issubset(D['phase_02'].keys())

    def test_to_dict_interleaved_categorized_by_data_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['data', 'element'])

        b = set(['net_01', 'phase_01', 'phase_02'])
        assert set(D.keys()) == b

        d = set(['pore', 'throat'])
        assert d.issubset(D['net_01']['labels'].keys())
        assert d.issubset(D['net_01']['properties'].keys())
        assert d.issubset(D['phase_01']['labels'].keys())
        assert d.issubset(D['phase_01']['properties'].keys())
        assert d.issubset(D['phase_02']['labels'].keys())
        assert d.issubset(D['phase_02']['properties'].keys())

    def test_to_dict_categorize_by_project(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['project'])
        assert 'proj_01' in D.keys()

    def test_from_dict_interleaved_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=True,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 2
        assert len(proj.geometries().values()) == 0
        assert len(proj.phases().values()) == 1
        assert len(proj.physics().values()) == 0

    def test_from_dict_interleaved_not_categorized(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=True,
                         categorize_by=[])
        proj = Dict.from_dict(D)
        assert len(proj) == 2
        assert len(proj.geometries().values()) == 0
        assert len(proj.phases().values()) == 0
        assert len(proj.physics().values()) == 0

    def test_from_dict_not_interleaved_flatted_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=True, interleave=False,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 6
        assert len(proj.geometries().values()) == 2
        assert len(proj.phases().values()) == 1
        assert len(proj.physics().values()) == 2

    def test_from_dict_not_interleaved_not_flatted_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 6
        assert len(proj.geometries().values()) == 2
        assert len(proj.phases().values()) == 1
        assert len(proj.physics().values()) == 2

    def test_from_dict_not_interleaved_not_flatted_cat_by_obj_data_elem(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=['object', 'element', 'data'])
        # Ensure that data and element categorizations are ripped out
        proj = Dict.from_dict(D)
        assert len(proj) == 6
        assert len(proj.geometries().values()) == 2
        assert len(proj.phases().values()) == 1
        assert len(proj.physics().values()) == 2

    def test_from_dict_not_interleaved_not_flatted_not_categorized(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=[])
        proj = Dict.from_dict(D)
        assert len(proj) == 6
        assert len(proj.geometries().values()) == 0
        assert len(proj.phases().values()) == 0
        assert len(proj.physics().values()) == 0

    def test_save_and_load(self, tmpdir):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=[])
        fname = tmpdir.join('test.dct')
        Dict.export_data(dct=D, filename=fname)
        dct = Dict.import_data(filename=fname)
        assert len(dct.keys()) == 2
        os.remove(fname)


if __name__ == '__main__':

    t = DictTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
