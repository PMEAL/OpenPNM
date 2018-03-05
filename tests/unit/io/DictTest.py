import openpnm as op
from openpnm.io import Dict
import scipy as sp
import pytest


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
        Ts = self.net.find_neighbor_throats(pores=Ps, mode='intersection')
        self.geo_2 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_2['pore.boo'] = 1
        self.geo_2['throat.boo'] = 1

        self.phase_1 = op.phases.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phases.GenericPhase(network=self.net)
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
        ws = op.core.Workspace()
        ws.clear()

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

        a = set(['net_01', 'geo_01', 'geo_02', 'phase_01', 'phys_01',
                 'phys_02', 'phase_02', 'phys_03', 'phys_04'])

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

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

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

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

        # Ensure categorized by object
        assert a == set(D.keys())

        # Ensure flatten, which occurs when categorized by object
        assert set(['geo_01', 'geo_02']).issubset(D['geometry'].keys())
        assert set(['phase_01', 'phase_02']).issubset(D['phase'].keys())
        assert set(['phys_01', 'phys_02']).issubset(D['physics'].keys())

    def test_to_dict_not_flat_not_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['data'])

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

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

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
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

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

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

        # Ensure categorized by element
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

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

        # Check if categorized by object
        assert a == set(D.keys())

        # Ensure it's flattened, which occurs when categorized by object
        assert set(['geo_01', 'geo_02']).isdisjoint(D['network']['net_01'].keys())
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase']['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase']['phase_02'].keys())

        # Check if categorized by data
        assert c.issubset(D['network']['net_01'].keys())
        assert c.issubset(D['phase']['phase_01'].keys())
        assert c.issubset(D['phase']['phase_02'].keys())
        assert c.issubset(D['geometry']['geo_01'].keys())
        assert c.issubset(D['geometry']['geo_02'].keys())
        assert c.issubset(D['physics']['phys_01'].keys())
        assert c.issubset(D['physics']['phys_02'].keys())
        assert c.issubset(D['physics']['phys_03'].keys())
        assert c.issubset(D['physics']['phys_04'].keys())

        # Ensure it's categorized by element
        assert d.issubset(D['network']['net_01']['labels'].keys())
        assert d.issubset(D['phase']['phase_01']['properties'].keys())
        assert d.issubset(D['phase']['phase_01']['labels'].keys())
        assert d.issubset(D['phase']['phase_02']['properties'].keys())
        assert d.issubset(D['phase']['phase_02']['labels'].keys())
        assert d.issubset(D['geometry']['geo_01']['properties'].keys())
        assert d.issubset(D['geometry']['geo_01']['labels'].keys())
        assert d.issubset(D['geometry']['geo_02']['properties'].keys())
        assert d.issubset(D['geometry']['geo_02']['labels'].keys())
        assert d.issubset(D['physics']['phys_01']['properties'].keys())
        assert d.issubset(D['physics']['phys_01']['labels'].keys())
        assert d.issubset(D['physics']['phys_02']['properties'].keys())
        assert d.issubset(D['physics']['phys_02']['labels'].keys())
        assert d.issubset(D['physics']['phys_03']['properties'].keys())
        assert d.issubset(D['physics']['phys_03']['labels'].keys())
        assert d.issubset(D['physics']['phys_04']['properties'].keys())
        assert d.issubset(D['physics']['phys_04']['labels'].keys())

    def test_to_dict_not_flat_not_interleaved_cat_by_element_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=False,
                         categorize_by=['element', 'object'])

        a = set(['network', 'phase', 'physics', 'geometry'])
        b = set(['net_01', 'phase_01', 'phase_02'])
        c = set(['labels', 'properties'])
        d = set(['pore', 'throat'])

        # Check if categorized by object
        assert a == set(D.keys())

        # Ensure it's flattened, which occurs when categorized by object
        assert set(['geo_01', 'geo_02']).isdisjoint(D['network']['net_01'].keys())
        assert set(['phys_01', 'phys_02']).isdisjoint(D['phase']['phase_01'].keys())
        assert set(['phys_03', 'phys_04']).isdisjoint(D['phase']['phase_02'].keys())

        # Ensure it's categorized by element
        assert d.issubset(D['network']['net_01'].keys())
        assert d.issubset(D['phase']['phase_01'].keys())
        assert d.issubset(D['phase']['phase_01'].keys())
        assert d.issubset(D['phase']['phase_02'].keys())
        assert d.issubset(D['phase']['phase_02'].keys())
        assert d.issubset(D['geometry']['geo_01'].keys())
        assert d.issubset(D['geometry']['geo_01'].keys())
        assert d.issubset(D['geometry']['geo_02'].keys())
        assert d.issubset(D['geometry']['geo_02'].keys())
        assert d.issubset(D['physics']['phys_01'].keys())
        assert d.issubset(D['physics']['phys_01'].keys())
        assert d.issubset(D['physics']['phys_02'].keys())
        assert d.issubset(D['physics']['phys_02'].keys())
        assert d.issubset(D['physics']['phys_03'].keys())
        assert d.issubset(D['physics']['phys_03'].keys())
        assert d.issubset(D['physics']['phys_04'].keys())
        assert d.issubset(D['physics']['phys_04'].keys())

    def test_to_dict_flat_not_interleaved_categorized_by_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['element'])

    def test_to_dict_flat_not_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['data'])

    def test_to_dict_flat_not_interleaved_categorized_by_data_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=True, interleave=False,
                         categorize_by=['data', 'element'])

    def test_to_dict_interleaved_categorized_by_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['element'])

    def test_to_dict_interleaved_categorized_by_data(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['data'])

    def test_to_dict_interleaved_categorized_by_data_element(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1, self.phase_2],
                         flatten=False, interleave=True,
                         categorize_by=['data', 'element'])

    def test_from_dict_interleaved_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=True,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 2

    def test_from_dict_not_interleaved_flatted_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=True, interleave=False,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 6

    def test_from_dict_not_interleaved_not_flatted_categorized_by_object(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=['object'])
        proj = Dict.from_dict(D)
        assert len(proj) == 6

    def test_from_dict_not_interleaved_not_flatted_not_categorized(self):
        D = Dict.to_dict(network=self.net, phases=[self.phase_1],
                         flatten=False, interleave=False,
                         categorize_by=[])
        proj = Dict.from_dict(D)
        assert len(proj) == 6

if __name__ == '__main__':

    t = DictTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
