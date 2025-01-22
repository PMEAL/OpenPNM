import pnmlib as pl
import numpy as np
import pytest


class GenericTest:

    hr = '―' * 78

    def __init__(self):
        print(self.hr)

    def setup_class(self):
        self.prj = {}
        pn = pl.generators.cubic([3, 3, 1])
        self.prj['network'] = pn

    def run_all(self):
        self.setup_class()
        for item in self.__dir__():
            if item.startswith('test'):
                print(f"Running test: {item}")
                self.__getattribute__(item)()


class CoreTests(GenericTest):

    def test_count(self):
        assert pl.core.count(self.prj['network'], element='pore') == 9
        assert pl.core.count(self.prj['network'], element='throat') == 12


# pl.inspect.tree(d, hide=None)


if __name__ == '__main__':

    t = CoreTests()
    t.run_all()
