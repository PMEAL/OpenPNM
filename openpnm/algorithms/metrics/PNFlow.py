import subprocess
from pathlib import Path


class PNFlow:
    r"""
    run pnflow.exe file on a network

    Notes
    -----
    This function only runs on Windows since the Windows compatible binary is
    provided by the Imperial College team.
    """
    @classmethod
    def run(cls, prefix, path='./'):
        r"""
        run pnflow.exe on a network

        Parameters
        ----------
        prefix : string
            The prefix of the filenames (i.e. 'prefix_node1.dat')

        Notes
        -----
        All the required files including prefix_link1, prefix_link2,
        prefix_node1, prefix_node2, prefix.dat and pnflow.exe file
        should be in a same directory as the exectuable.
        """

        # Create dat file
        s = f"TITLE  {prefix}; \n" \
            f"writeStatistics true; \n" \
            f"NETWORK  F {prefix} \n" \
            f"cycle1:   0.00     2.00E+05     0.05      T     T; \n" \
            f"cycle2:   1.00    -2.00E+05     0.05      T     T; \n" \
            f"cycle3:   0.00     2.00E+05     0.05      T     T; \n" \
            f"cycle1_BC:   T     F       T       T      DP    1.00  2.00; \n" \
            f"cycle2_BC:   T     F       T       T      DP    2.00  1.00; \n" \
            f"cycle3_BC:   T     F       T       T      DP    1.00  2.00; \n" \
            f"CALC_BOX:  0.1 0.9; \n" \
            f"INIT_CONT_ANG:   1   0   10  -0.2    -3.0   rand   0.0; \n" \
            f"EQUIL_CON_ANG:   4   30   50   -0.2   -3.0   rand   0.0; \n" \
            f"Water  0.001               1.2                1000.0; \n" \
            f"Oil   0.001               1000.0              1000.0; \n" \
            f"ClayResistivity            2.0 ; \n" \
            f"WaterOilInterface          0.03 ; \n" \
            f"DRAIN_SINGLETS: T; \n" \
            f"SAT_CONTROL: F; \n"

        file = Path(path + prefix + '.dat')
        with open(file, 'w') as f:
            f.write(s)

        exe = Path(path + 'pnflow.exe')
        subprocess.Popen([exe, file])
