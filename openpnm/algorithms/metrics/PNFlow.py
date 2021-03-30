import subprocess


class PNFlow:
    r"""
    run pnflow.exe file on a network

    Notes
    -----
    This function only runs on Windows since the Windows compatible binary is
    provided by the Imperial College team.
    """
    @classmethod
    def run(cls, prefix, path_to_exe):
        r"""
        run pnflow.exe file on a network
        ----------
        prefix : string
            The prefix of the filenames (i.e. 'prefix_node1.dat')

        path_to_exe : string
            Path to the pnflow .exe file. See Notes

        Notes
        -----
        all the required files including prefix_link1, prefix_link2,
        prefix_node1, prefix_node2, prefix.dat and pnflow.exe file
        should be in a same directory
        """
        subprocess.Popen([path_to_exe, prefix+'.dat'])
