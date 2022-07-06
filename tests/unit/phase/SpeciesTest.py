import openpnm as op


class SpeciesTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = SpeciesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
