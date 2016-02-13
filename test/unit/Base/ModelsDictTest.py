import pytest
import OpenPNM
import scipy as sp
from OpenPNM.Base import ModelsDict
from OpenPNM.Base.__ModelsDict__ import ModelWrapper


class ModelsDictTest:

    class ModelWrapperTest:
        def test_str_with_none_model(self):
            actual_string = ModelWrapper(model=None).__str__()
            expected_string = 'No model specified.'
            assert expected_string == actual_string

        def test_find_master_with_more_than_one_master(self):
            pn = OpenPNM.Network.Cubic(shape=[3, 3, 3])
            geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
                                                    throats=pn.Ts)
            geom.models.add(propname='pore.seed',
                            model=OpenPNM.Geometry.models.pore_misc.random)
            pn2 = OpenPNM.Network.Cubic(shape=[3, 3, 3])
            geom2 = OpenPNM.Geometry.GenericGeometry(network=pn2, pores=pn.Ps,
                                                     throats=pn.Ts)
            geom2.models = geom.models
            with pytest.raises(Exception):
                geom2.models['pore.seed'].run()

    class ModelsDictTest:
        def setup_class(self):
            self.models_dict = ModelsDict()

        def test_str(self):
            actual_string = self.models_dict.__str__()
            expected_string = \
                '------------------------------------------------------------\n' + \
                '#     Property Name                  Regeneration Mode\n' + \
                '------------------------------------------------------------\n' + \
                '------------------------------------------------------------'
            assert expected_string == actual_string

        def test_keys(self):
            assert self.models_dict.keys() == []
            self.models_dict['test'] = {'model': 'test content'}
            assert self.models_dict.keys() == ['test']

        def test_add_with_no_master(self):
            with pytest.raises(Exception):
                self.models_dict.add('test model', {})

        def test_copy_modelsdict_to_new_object(self):
            pn = OpenPNM.Network.Cubic(shape=[3, 3, 3])
            geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
                                                    throats=pn.Ts)
            geom.models.add(propname='pore.seed',
                            model=OpenPNM.Geometry.models.pore_misc.random)
            pn2 = OpenPNM.Network.Cubic(shape=[3, 3, 3])
            geom2 = OpenPNM.Geometry.GenericGeometry(network=pn2, pores=pn.Ps,
                                                     throats=pn.Ts)
            geom2.models = geom.models
            with pytest.raises(Exception):
                geom2.regenerate()
            geom2.models = geom.models.copy()
            geom2.regenerate()
            assert 'pore.seed' in geom2

        def test_find_master_with_more_than_one_master(self):
            pass

        def test_reorder(self):
            pn = OpenPNM.Network.Cubic(shape=[5, 5, 5])
            geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
                                                    throats=pn.Ts)
            geom.models.add(propname='throat.blah',
                            model=OpenPNM.Geometry.models.throat_misc.random,
                            regen_mode='deferred')
            geom.models.add(propname='throat.seed',
                            model=OpenPNM.Geometry.models.throat_misc.neighbor,
                            pore_prop='pore.seed',
                            mode='min',
                            regen_mode='deferred')
            geom.models.add(propname='pore.seed',
                            model=OpenPNM.Geometry.models.pore_misc.random,
                            seed=None,
                            regen_mode='deferred')
            with pytest.raises(Exception):
                geom.regenerate()
            geom.models.reorder({'pore.seed': 1, 'throat.seed': 2})
            geom.regenerate()
            assert 'throat.seed' in geom
            geom.models.move_to_end('throat.blah')
            string = geom.models.__str__()
            expected = \
                '------------------------------------------------------------\n' + \
                '#     Property Name                  Regeneration Mode\n' + \
                '------------------------------------------------------------\n' + \
                '1     pore.seed                      deferred            \n' + \
                '2     throat.seed                    deferred            \n' + \
                '3     throat.blah                    deferred            \n' + \
                '------------------------------------------------------------'
            assert string == expected

        def test_remove(self):
            pn = OpenPNM.Network.Cubic(shape=[5, 5, 5])
            geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
                                                    throats=pn.Ts)
            geom.models.add(propname='pore.seed',
                            model=OpenPNM.Geometry.models.pore_misc.random,
                            seed=None,
                            regen_mode='constant')
            assert 'pore.seed' in geom
            assert 'pore.seed' in geom.models
            geom.remove('pore.seed')
            assert 'pore.seed' not in geom
            assert 'pore.seed' not in geom.models

        def test_changing_regen_mode(self):
            pn = OpenPNM.Network.Cubic(shape=[5, 5, 5])
            geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps,
                                                    throats=pn.Ts)
            geom.models.add(propname='pore.seed',
                            model=OpenPNM.Geometry.models.pore_misc.random,
                            seed=None,
                            regen_mode='constant')
            a = sp.copy(geom['pore.seed'])
            geom.regenerate()
            assert sp.all(a == geom['pore.seed'])
            geom.models['pore.seed']['regen_mode'] = 'normal'
            geom.regenerate()
            assert not sp.all(a == geom['pore.seed'])
