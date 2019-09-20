import os
import subprocess
import logging
import openpnm as op

rootdir = os.path.split(os.path.split(op.__file__)[0])[0]
examples_dir = os.path.join(rootdir, 'examples')
test_dir = os.path.join(examples_dir, 'topology')


class TopoTest():

    def setup_class(self):
        pass

    def _run_shell_command(self, command_line_args):
        try:
            proc = subprocess.run(command_line_args, timeout=600)
        except (OSError, subprocess.CalledProcessError) as exception:
            logging.info('Exception occured: ' + str(exception))
            logging.info('Subprocess failed')
            return False
        else:
            # no exception was raised
            logging.info('Subprocess finished')
        return proc.returncode == 0

    def _notebook_run(self, path):
        """Execute a notebook via nbconvert and collect output.
           :returns (parsed nb object, execution errors)
        """
        dirname, __ = os.path.split(path)
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=360",
                "--output", "temp_output.ipynb", path]
        rc = self._run_shell_command(args)
        print(path, rc)
        print('-'*30)
        if rc:
            os.remove(os.path.join(dirname, "temp_output.ipynb"))
        return rc

    def test_Generate_Dual_Cubic_Lattice(self):
        nbook = os.path.join(test_dir, 'Generate Dual Cubic Lattice.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Stitching_Networks_Together(self):
        nbook = os.path.join(test_dir, 'Stitching Networks Together.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Using_Queries_Level_1(self):
        nbook = os.path.join(test_dir, 'Using Queries - Level 1.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Various_Cubic_Networks(self):
        nbook = os.path.join(test_dir, 'Various Cubic Networks.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Voronoi_Fibers(self):
        nbook = os.path.join(test_dir, 'Voronoi Fibers.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Topology_Extension(self):
        nbook = os.path.join(test_dir, 'Topology Extension.ipynb')
        print(nbook)
        rc = self._notebook_run(nbook)
        assert rc

if __name__ == '__main__':
    t = TopoTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
