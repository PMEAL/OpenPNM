import openpnm as op
import numpy as np
import openpnm.models as mods
import pytest


class ModelsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_models_dict_print(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.StickAndBall(network=net, pores=net.Ps,
                                       throats=net.Ts)
        s = geo.models.__str__().split('\n')
        assert len(s) == 69
        assert s.count('―'*85) == 15

    def test_regenerate_models(self):
        a = len(self.geo.props())
        assert a == 16
        self.geo.clear(mode='props')
        a = len(self.geo.props())
        assert a == 0
        self.geo.regenerate_models()
        a = len(self.geo.props())
        assert a == 16

    def test_dependency_graph(self):
        phase = op.phases.GenericPhase(network=self.net)
        phase.add_model(propname="pore.foo", model=op.models.misc.constant, value=1.0)
        phys = op.physics.GenericPhysics(network=self.net,
                                         phase=phase,
                                         geometry=self.geo)
        phys.add_model(propname="pore.baz", model=op.models.misc.constant, value=0.0)

        def mymodel(target, foo="pore.foo", baz="pore.baz"):
            return 0.0
        phys.add_model(propname="pore.bar_depends_on_foo_and_baz", model=mymodel)
        dg = phys.models.dependency_graph()
        assert ["pore.baz", "pore.bar_depends_on_foo_and_baz"] in dg.edges()
        assert ["pore.foo", "pore.bar_depends_on_foo_and_baz"] not in dg.edges
        dg = phys.models.dependency_graph(deep=True)
        assert ["pore.baz", "pore.bar_depends_on_foo_and_baz"] in dg.edges
        assert ["pore.foo", "pore.bar_depends_on_foo_and_baz"] in dg.edges

    def test_dependency_list(self):
        prj = self.net.project
        prj.purge_object(self.geo)
        geom = op.geometry.GenericGeometry(network=self.net,
                                           pores=self.net.Ps)

        geom.add_model(propname='pore.volume',
                       model=mods.geometry.pore_volume.sphere,
                       pore_diameter='pore.diameter',
                       regen_mode='deferred')

        geom.add_model(propname='pore.diameter',
                       model=mods.misc.product,
                       prop1='pore.max_size',
                       prop2='pore.seed',
                       regen_mode='deferred')

        geom.add_model(propname='pore.area',
                       model=mods.geometry.pore_cross_sectional_area.sphere,
                       pore_diameter='pore.diameter',
                       regen_mode='deferred')

        geom.add_model(propname='pore.seed',
                       model=mods.misc.random,
                       element='pore',
                       num_range=[0, 0.1],
                       seed=None)

        tree = np.asarray(geom.models.dependency_list())
        pos_v = np.argwhere(tree == 'pore.volume').flatten()[0]
        pos_d = np.argwhere(tree == 'pore.diameter').flatten()[0]
        pos_a = np.argwhere(tree == 'pore.area').flatten()[0]
        pos_s = np.argwhere(tree == 'pore.seed').flatten()[0]
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

        with pytest.raises(Exception):
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

    def test_regenerate_models_on_phase_with_deep(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
        phase = op.phases.Water(network=pn)
        phys = op.physics.Standard(network=pn, phase=phase, geometry=geo)
        phase.clear(mode='model_data')
        phys.clear()
        assert len(phys) == 2  # Only pore and throat.all remain
        phase.regenerate_models(propnames=None, deep=False)
        assert len(phys) == 2  # Still only pore and throat.all
        phase.regenerate_models(propnames=None, deep=True)
        assert len(phys) > 2  # Phys models are regenerated by phase regen

    def test_regenerate_models_on_physics_with_deep(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
        phase = op.phases.Water(network=pn)
        phys = op.physics.Standard(network=pn, phase=phase, geometry=geo)
        len_phase = 23
        phase.clear(mode='model_data')
        phys.clear()
        ws = op.Workspace()
        loglevel = ws.settings["loglevel"]
        ws.settings["loglevel"] = 50
        assert len(phys) == 2
        assert len(phase) == 13
        phys.regenerate_models(propnames=None, deep=False)
        assert len(phys) == 14
        # Note that 2 new models were added to the phase during interpolation
        assert len(phase) < len_phase
        phase.clear(mode='model_data')
        assert len(phase) == 13
        phys.regenerate_models(propnames=None, deep=True)
        assert len(phase) < len_phase
        ws.settings["loglevel"] = loglevel

    def test_regenerate_models_on_network_with_deep(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
        a = len(pn.props())
        pn.clear()
        pn.regenerate_models()
        assert len(pn.props()) == a  # pn has NO models
        b = len(geo.props())
        geo.clear(mode='model_data')
        pn.regenerate_models(deep=False)
        assert len(geo.props()) == 0
        pn.regenerate_models(deep=True)
        assert len(geo.props()) == b

    def test_regen_mode_default_value(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=1e-4)
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts,
                                        settings={'regen_mode': 'deferred'})
        assert len(geo.props()) == 0
        geo.regenerate_models()
        assert len(geo.props()) == 16

    def test_automatic_running_on_models_when_missing_data(self):
        pn = op.network.Cubic(shape=[3, 3, 3], spacing=1e-4)
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts,
                                        settings={'regen_mode': 'deferred'})
        assert len(geo) == 2
        _ = geo['pore.seed']
        assert len(geo) == 3


if __name__ == '__main__':

    t = ModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
