import os
import subprocess
import logging
import openpnm as op

rootdir = os.path.split(os.path.split(op.__file__)[0])[0]
examples_dir = os.path.join(rootdir, 'examples')
test_dir = os.path.join(examples_dir, 'io_and_visualization')


class IOVisTest():

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

    def test_Quick_Plotting_in_OpenPNM(self):
        nbook = os.path.join(test_dir, 'Quick Plotting in OpenPNM.ipynb')
        rc = self._notebook_run(nbook)
        assert rc

    def test_Statoil_Import_and_Permeability_Calculation(self):
        fname = 'Statoil Import and Permeability Calculation.ipynb'
        nbook = os.path.join(test_dir, fname)
        rc = self._notebook_run(nbook)
        assert rc


if __name__ == '__main__':
    t = IOVisTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
