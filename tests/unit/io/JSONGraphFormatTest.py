import openpnm as op
import py


class JSONGraphFormatTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = JSONGraphFormatTest()
    self = t  # For interacting with the tests at the command line
    tmpdir = py.path.local()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=tmpdir)
