import openpnm as op
import numpy as np


class PandasTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2], name='bob')
        self.net['pore.boo'] = 1
        self.net['throat.boo'] = 1
        self.phase_1 = op.phase.Phase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phase.Phase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2
        self.phase_1['pore.baz'] = 11
        self.phase_1['throat.baz'] = 11
        self.phase_2['pore.baz'] = 12
        self.phase_2['throat.baz'] = 12

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_project_to_dataframe_not_joined(self):
        df = op.io.project_to_pandas(project=self.net.project,
                                     join=False)
        assert len(df['pore'].keys()) == 21
        assert len(df['throat'].keys()) == 10

    def test_project_to_dataframe_joined(self):
        df = op.io.project_to_pandas(project=self.net.project,
                                     join=True)
        assert len(df.keys()) == 31
        assert np.isnan(df['bob.pore.coords[0]']).sum() > 0

    def test_network_to_dataframe(self):
        df = op.io.network_to_pandas(network=self.net)
        assert len(df.keys()) == 2
        assert len(df['pore'].keys()) == 11
        assert len(df['throat'].keys()) == 4


if __name__ == '__main__':
    t = PandasTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
