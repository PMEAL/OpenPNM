import numpy as np
import openpnm as op


class GeneratorTest:

    def setup_class(self):
        self.ws = op.Workspace()




if __name__ == '__main__':

    t = GeneratorTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
