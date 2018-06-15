import openpnm as op
import py
import os
import pytest
from pathlib import Path


class GenericIOTest:

    def setup_class(self):
        self._count = 0
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase2 = op.phases.GenericPhase(network=self.net)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()

    def count(self):
        self._count += 1
        return self._count

    @staticmethod
    def create_write_and_close_file(fname):
        r"""
        Putting the following code into it's own function saves a lot of
        repitition, but means the assert errors will not come directly from
        the offending test.  The pytest output gives a traceback which leads
        to the correct spot, so this will suffice.
        """
        # Make sure file does not exist yet
        assert not fname.is_file()
        # Make file and write a line to it
        with open(fname, 'w') as f:
            f.write('write a line in file')
        # Confirm file is present
        assert fname.is_file()
        # Remove file for good measure
        os.remove(fname)

    def test_parse_filename_dot_extension_no_path(self, tmpdir):
        filename = 'test'+str(self.count())+'.ext'
        fname = op.io.GenericIO._parse_filename(filename=filename)
        self.create_write_and_close_file(fname)

    def test_parse_filename_arg_extension_no_path(self):
        filename = 'test'+str(self.count())
        fname = op.io.GenericIO._parse_filename(filename=filename, ext='ext')
        self.create_write_and_close_file(fname)

    def test_parse_filename_arg_and_dot_extension_no_path(self):
        filename = 'test'+str(self.count())+'.ext'
        fname = op.io.GenericIO._parse_filename(filename=filename, ext='other')
        self.create_write_and_close_file(fname)

    def test_parse_filename_dot_extension_with_path(self, tmpdir):
        filename = 'test'+str(self.count())+'.ext'
        filename = os.path.join(tmpdir.dirname, filename)
        fname = op.io.GenericIO._parse_filename(filename=filename)
        self.create_write_and_close_file(fname)

    def test_parse_filename_path_in_filename_dot_extension(self, tmpdir):
        filename = Path(tmpdir.join('test'+str(self.count())+'.ext'))
        fname = op.io.GenericIO._parse_filename(filename=filename)
        self.create_write_and_close_file(fname)

    def test_parse_args_no_input_lists(self):
        proj, net, phases = op.io.GenericIO._parse_args(network=self.net,
                                                        phases=self.phase)
        assert isinstance(proj, list)
        assert isinstance(net, list)
        assert isinstance(phases, list)

    def test_parse_args_all_input_lists(self):
        proj, net, phases = op.io.GenericIO._parse_args(network=[self.net],
                                                        phases=[self.phase])
        assert isinstance(proj, list)
        assert isinstance(net, list)
        assert isinstance(phases, list)

    def test_parse_args_repeated(self):
        proj1, net1, phases1 = op.io.GenericIO._parse_args(network=[self.net],
                                                           phases=[self.phase])
        proj2, net2, phases2 = op.io.GenericIO._parse_args(network=[net1],
                                                           phases=[phases1])
        assert proj1 == proj2
        assert net1 == net2
        assert phases1 == phases2

    def test_save(self):
        with pytest.raises(NotImplementedError):
            op.io.GenericIO.save()

    def test_load(self):
        with pytest.raises(NotImplementedError):
            op.io.GenericIO.load()


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = GenericIOTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
