import openpnm as op
import numpy as np
import pytest


class HealthCheckTest:

    def setup_class(self):
        self.ws = op.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[2, 2, 2])

    def check_data_health(self):
        self.net.update({'pore.test': np.array([1, 2, 3, 4, 5, 6])})
        a = op.utils.check_data_health(self.net)
        assert a.health == False
        assert a['pore.test'] != []
        assert a['pore.coords'] == []
        assert a['throat.conns'] == []
        self.net['pore.test'] = 1.0
        self.net['pore.test'][0] = np.nan
        a = op.utils.check_data_health(self.net)
        assert a['pore.test'] != []


if __name__ == '__main__':

    t = HealthCheckTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
