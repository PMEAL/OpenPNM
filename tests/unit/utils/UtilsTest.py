import openpnm as op
import scipy as sp
import pytest


class UtilsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_tic_toc(self):
        with pytest.raises(Exception):
            op.utils.toc()
        op.utils.tic()
        t1 = op.utils.toc()
        assert t1 is None
        t2 = op.utils.toc(quiet=True)
        assert t2 >= 0

    def test_nested_dict(self):
        d = op.utils.NestedDict()
        d['top']['middle']['bottom'] = 1
        assert d == {'top': {'middle': {'bottom': 1}}}
        s = d.__str__()
        assert s == '-top\n--middle\n---bottom\n'
        a = d.to_dict()
        assert type(a) is dict

    def test_printable_list(self):
        L = op.utils.PrintableList(['item1', 'item2', 'item2'])
        s = L.__str__().split('\n')
        assert len(s) == 5

    def test_printable_dict(self):
        D = op.utils.PrintableDict(**{'item1': 1, 'item2': 2, 'item3': sp.array([1, 2])})
        s = D.__str__().split('\n')
        assert len(s) == 7
        r = D.__repr__()
        assert r == "{'item1': 1, 'item2': 2, 'item3': array([1, 2])}"


if __name__ == '__main__':

    t = UtilsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
