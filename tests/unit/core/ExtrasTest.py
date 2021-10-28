import importlib
import openpnm as op
from openpnm import Workspace
from types import ModuleType
ws = Workspace()


class ExtrasTest:

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_initialize_GenericNetwork_without_args(self):
        net = op.network.GenericNetwork()
        assert set(net.keys()) == set(['pore.all', 'throat.all'])
        assert net.Np == 0
        assert net.Nt == 0

    def test_initialize_GenericGeometry_without_args(self):
        obj = op.geometry.GenericGeometry()
        assert set(obj.keys()) == set(['pore.all', 'throat.all'])
        assert obj.Np == 0
        assert obj.Nt == 0
        assert len(obj.project) == 1

    def test_initialize_StickAndBall_without_args(self):
        obj = op.geometry.SpheresAndCylinders(settings={'freeze_models': True})
        assert set(obj.keys()) == set(['pore.all', 'throat.all'])
        assert len(obj.models.keys()) > 0
        assert obj.Np == 0
        assert obj.Nt == 0
        assert len(obj.project) == 1

    def test_initialize_GenericPhase_without_args(self):
        obj = op.phases.GenericPhase()
        assert obj.Np == 0
        assert obj.Nt == 0
        assert len(obj.project) == 1

    def test_initialize_Air_without_args(self):
        obj = op.phases.Air(settings={'freeze_models': True})
        assert len(obj.keys()) > 4
        assert len(obj.models.keys()) > 0
        assert obj.Np == 0
        assert obj.Nt == 0
        assert len(obj.project) == 1

    def test_initialize_GenericPhysics_without_args(self):
        obj = op.physics.GenericPhysics()
        assert set(obj.keys()) == set(['pore.all', 'throat.all'])
        assert obj.Np == 0
        assert obj.Nt == 0
        assert len(obj.project) == 1

    def test_init_Standard_physics_without_args(self):
        obj = op.physics.Standard(settings={'freeze_models': True})
        assert len(obj.models) > 0

    def test_init_geometris_without_args(self):
        obj = op.geometry.SpheresAndCylinders(settings={'freeze_models': True})
        assert len(obj.models) > 0

    def test_init_phases_without_args(self):
        phases = [p for p in dir(op.phases) if not p.startswith('__')]
        mod = importlib.import_module('openpnm.phases')
        for p in phases:
            clss = getattr(mod, p)
            if not isinstance(clss, ModuleType):
                clss(settings={'freeze_models': True})

    def test_init_geometry_with_either_network_or_project(self):
        ws.clear()
        ws.settings['loglevel'] = 50
        classes = [c for c in dir(op.geometry) if not c.startswith('__')]
        mod = importlib.import_module('openpnm.geometry')
        for c in classes:
            clss = getattr(mod, c)
            if not isinstance(clss, ModuleType):
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(project=net.project)
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(network=net)

    def test_init_phase_with_either_network_or_project(self):
        ws.clear()
        ws.settings['loglevel'] = 50
        classes = [c for c in dir(op.phases) if not c.startswith('__')]
        mod = importlib.import_module('openpnm.phases')
        for c in classes:
            clss = getattr(mod, c)
            if not isinstance(clss, ModuleType):
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(project=net.project)
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(network=net)

    def test_init_physics_with_either_network_or_project(self):
        ws.clear()
        ws.settings['loglevel'] = 50
        classes = [c for c in dir(op.physics) if not c.startswith('__')]
        mod = importlib.import_module('openpnm.physics')
        for c in classes:
            clss = getattr(mod, c)
            if not isinstance(clss, ModuleType):
                net = op.network.Cubic(shape=[2, 2, 2])
                # Eventually we want to create physics objects without a phase
                # but for now this is necessary.
                phase = op.phases.GenericPhase(network=net)
                clss(project=net.project, phase=phase)
                phase = op.phases.GenericPhase(network=net)
                clss(network=net, phase=phase)

    def test_init_algorithm_with_either_network_or_project(self):
        ws.clear()
        ws.settings['loglevel'] = 50
        classes = [c for c in dir(op.algorithms) if not c.startswith('__')]
        mod = importlib.import_module('openpnm.algorithms')
        # The following 3 should be fixed to accept no phase arguments
        classes.pop(classes.index('NernstPlanck'))
        classes.pop(classes.index('NernstPlanckMultiphysicsSolver'))
        classes.pop(classes.index('TransientNernstPlanckMultiphysicsSolver'))
        # metrics should be ignored
        classes.pop(classes.index('metrics'))
        for c in classes:
            clss = getattr(mod, c)
            if not isinstance(clss, ModuleType):
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(project=net.project)
                net = op.network.Cubic(shape=[2, 2, 2])
                clss(network=net)


if __name__ == '__main__':

    t = ExtrasTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
