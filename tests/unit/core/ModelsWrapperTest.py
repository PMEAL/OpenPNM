import pytest
import openpnm as op
import numpy as np


class ModelsTest:

    def setup_class(self):
        self.net = op.network.Demo([4, 4, 1])

    def test_call(self):
        a = self.net['pore.diameter']
        b = self.net.models['pore.diameter@all']()
        assert np.all(a == b)

    def test_find_target(self):
        a = self.net.models['pore.diameter@all'].target
        assert a is self.net

    def test_name(self):
        a = 'pore.diameter@all'
        assert a == self.net.models['pore.diameter@all'].name

    def test_propname(self):
        a = 'pore.diameter'
        assert a == self.net.models['pore.diameter@all'].propname

    def test_domain(self):
        a = 'pore.all'
        assert a == self.net.models['pore.diameter@all'].domain

    def test_run_model(self):
        a = self.net['pore.seed'].copy()
        b = self.net['pore.diameter'].copy()
        self.net.run_model('pore.seed')
        c = self.net['pore.seed'].copy()
        d = self.net['pore.diameter'].copy()
        assert not np.any(a == c)
        assert np.all(b == d)
        self.net.run_model('pore.diameter')
        e = self.net['pore.diameter'].copy()
        assert not np.any(b == e)




if __name__ == '__main__':

    t = ModelsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
