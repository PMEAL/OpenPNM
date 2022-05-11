import openpnm as op
import numpy as np


class PandasTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2], name='bob')
        self.net['pore.boo'] = 1
        self.net['throat.boo'] = 1
        self.phase_1 = op.phase.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phase.GenericPhase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2
        self.phase_1['pore.baz'] = 11
        self.phase_1['throat.baz'] = 11
        self.phase_2['pore.baz'] = 12
        self.phase_2['throat.baz'] = 12

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_to_dataframe_not_joined(self):
        df = op.io.to_pandas(network=self.net, phases=[self.phase_1],
                             join=False)
        assert len(df['pore'].keys()) == 16
        assert len(df['throat'].keys()) == 7

    def test_to_dataframe_joined(self):
        df = op.io.to_pandas(network=self.net, phases=[self.phase_1],
                             join=True)
        assert len(df.keys()) == 23
        assert np.isnan(df['bob.pore.coords[0]']).sum() > 0


if __name__ == '__main__':
    t = PandasTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
