import openpnm as op
import numpy as np
import openpnm.models as mods
import pytest
from testfixtures import LogCapture, ShouldRaise


class ModelsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()

    def test_models_dict_print(self):
        s = self.geo.models.__str__().split('\n')
        assert len(s) == 56
        assert s.count('―'*78) == 13

    def test_regenerate_models(self):
        a = len(self.geo.props())
        assert a == 11
        self.geo.clear(mode='props')
        a = len(self.geo.props())
        assert a == 0
        self.geo.regenerate_models()
        a = len(self.geo.props())
        assert a == 11

    def test_dependency_list(self):
        prj=self.net.project
        prj.purge_object(self.geo)
        geom = op.geometry.GenericGeometry(network=self.net,
                                           pores=self.net.Ps)

        geom.add_model(propname='pore.volume',
                       model=mods.geometry.pore_volume.sphere,
                       pore_diameter='pore.diameter')

        geom.add_model(propname='pore.diameter',
                       model=mods.misc.product,
                       prop1='pore.max_size',
                       prop2='pore.seed')

        geom.add_model(propname='pore.area',
                       model=mods.geometry.pore_area.sphere,
                       pore_diameter='pore.diameter')

        geom.add_model(propname='pore.seed',
                       model=mods.misc.random,
                       element='pore',
                       num_range=[0, 0.1],
                       seed=None)
        tree = np.asarray(geom.models.dependency_list())
        pos_v = np.argwhere(tree=='pore.volume').flatten()[0]
        pos_d = np.argwhere(tree=='pore.diameter').flatten()[0]
        pos_a = np.argwhere(tree=='pore.area').flatten()[0]
        pos_s = np.argwhere(tree=='pore.seed').flatten()[0]
        assert pos_v > pos_d
        assert pos_d > pos_s
        assert pos_a > pos_d
        self.geo = geom

    def test_dependency_list_circular(self):
        pn = self.net

        def chicken(target, prop='pore.egg'):
            return np.ones(target.Np)

        def egg(target, prop='pore.chicken'):
            return np.ones(target.Np)

        pn.add_model(propname='pore.chicken', model=chicken)
        pn.add_model(propname='pore.egg', model=egg)

        with ShouldRaise(Exception('Cyclic dependency found: pore.egg ->' +
                                   ' pore.chicken -> pore.egg')):
            pn.models.dependency_list()
        pn.remove_model('pore.chicken')
        pn.remove_model('pore.egg')

    def test_dependency_list_tri_circular(self):
        pn = self.net

        def rock(target, prop='pore.scissors'):
            return np.ones(target.Np)

        def scissors(target, prop='pore.paper'):
            return np.ones(target.Np)

        def paper(target, prop='pore.rock'):
            return np.ones(target.Np)

        pn.add_model(propname='pore.paper', model=paper)
        pn.add_model(propname='pore.scissors', model=scissors)
        pn.add_model(propname='pore.rock', model=rock)

        with pytest.raises(Exception):
            pn.models.dependency_list()


if __name__ == '__main__':

    t = ModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
